function [pulse_model, pLik] = getPulseLikelihood(NewPulses, PulseModel, scaleType, smoothGate)
%[pulse_model,Lik_pulse] = Z_2_pulse_model(pulse_model,new_pulses)
%USAGE
% assumes both pulse model and new pulses have the same sampling rate
% input: NewPulses: new pulses, and PulseModel: pulse model
%return pulse model & std etc and Lik of individual pulses given the model

%fit_pulse_model estimates only the fundamental frequency model using data
%that best fits this model. It then decimates the model and best fit data
%to build the second harmonic models for likelihood testing

%e.g.
%     pulse_model = 
%     fhM: [1x183 double]; mean first harmonic
%     shM: [1x168 double]; mean second harmonic
%     S_fhM: [1x183 double]; std first harmonic
%     S_shM: [1x168 double]; std second harmonic
if ~isempty(NewPulses)
    fprintf('Fitting to pulse model ')
    if ~exist('scaleType','var'); scaleType = 1; end
    if ~exist('smoothGate','var'); smoothGate = 0; end
    
    % collecting model info
    fhM = PulseModel.fhM;
    shM = PulseModel.shM;
    lengthfhM = length(fhM);
    lengthshM = length(shM);
    S_fhM = double(PulseModel.S_fhM); % std
    S_shM = double(PulseModel.S_shM); % std

    % collecting new pulses info and storing them into cells
    if ~iscell(NewPulses)
        if size(NewPulses, 1) == 251; NewPulses = NewPulses'; end
        NewPulses = mat2cell(NewPulses', 251, ones(size(NewPulses, 1), 1));
        if smoothGate
            NewPulses = cellfun(@(x) smooth(x, 10), NewPulses, 'UniformOutput', false);
        end
    end

    %grab samples, center, and pad
    n_samples = length(NewPulses);
    max_length = max(cellfun(@length, NewPulses));
    total_length = 2*max([lengthfhM  lengthshM max_length]);
    Z = zeros(n_samples, total_length);
    % I will assume all pulses are centered at the max power 
    % (given the updating of wc that I do based on this)
    NewCenter = round(total_length/2);
    for n = 1:n_samples
        X = NewPulses{n};
        Xlength = round(length(X)/2)-1;
        if mod(size(X,2),2)
            Z(n,[NewCenter - Xlength:NewCenter + Xlength]) = X(:);
        else
            Z(n,[NewCenter - Xlength:NewCenter + Xlength + 1]) = X(:);
        end
    end
    clear X Xlength

    %pad models / align to model (f and s) and phase reversed model
    pModFields = [];
    if ~isempty(fhM)
        pModFields = [pModFields, {'fhM','S_fhM'}]; 
    end
    if ~isempty(shM)
        pModFields = [pModFields, {'shM','S_shM'}];
    end
    % allocating aligned - scaled - trimmed new pulses to 1st and 2nd harminic
    Z2fhM = []; Z2shM = [];
    for pMod_idx = 1:numel(pModFields)
        storeVect = zeros(1, total_length);
        eval(['Xlength = round(length(',pModFields{pMod_idx},')/2)-1;'])
        if eval(['mod(size(', pModFields{pMod_idx}, ', 2), 2)'])
            eval(['storeVect(1,[NewCenter - Xlength:NewCenter + Xlength]) = ', pModFields{pMod_idx}, ';'])
        else
            eval(['storeVect(1,[NewCenter - Xlength:NewCenter + Xlength + 1]) = ', pModFields{pMod_idx}, ';'])
        end
        eval([pModFields{pMod_idx}, ' = storeVect;',])
        clear storeVect
        if strcmp(pModFields{pMod_idx}, 'fhM')
            % first harmonic
            Z2fhM = alignpulses2model(Z, fhM);
            Z2fhM = scaleZ2M(Z2fhM, fhM, scaleType);
            % Generate phase reversed model
            fhRM = -fhM;
            Z2fhRM = alignpulses2model(Z, fhRM);
            Z2fhRM = scaleZ2M(Z2fhRM, fhRM, scaleType);
        elseif strcmp(pModFields{pMod_idx}, 'shM')
            % second harmonic
            Z2shM = alignpulses2model(Z, shM);
            Z2shM = scaleZ2M(Z2shM, shM, scaleType);
            % Generate phase reversed model
            shRM = -shM;
            Z2shRM = alignpulses2model(Z, shRM);
            Z2shRM = scaleZ2M(Z2shRM, shRM, scaleType);
        end
    end
    clear Xlength pModFields pModFields

    % calculating chi-square
    chisq = nan(n_samples, 1);
    for n = 1:n_samples
        if ~isempty(fhM)
            %calc chi-square for first model
            chisq(n,1) = ...
                mean((Z2fhM(n, :) - fhM).^2./var(Z2fhM(n, :)));
            %calc chi-square for reversed model
            chisq(n,3) = ...
                mean((Z2fhRM(n, :) - fhRM).^2./var(Z2fhRM(n, :)));
        end
        if ~isempty(shM)
            %calc chi-square for second harmonic model
            chisq(n,2) = ...
                mean((Z2shM(n, :) - shM).^2./var(Z2shM(n, :)));
            %calc chi-square for reversed second harmonic model
            chisq(n,4) = ...
                mean((Z2shRM(n, :) - shRM).^2./var(Z2shRM(n, :)));
        end
    end
    [~, best_chisqr_idx] = nanmin(chisq, [], 2);

    %flip data that fits a reversed model better (columns 3 or 4)
    for n = 1:n_samples
        if best_chisqr_idx(n) > 2
            Z(n, :) = -Z(n, :);
        end
    end
    clear Z2fhRM Z2shRM  fhRM shRM n best_chisqr_idx chisq

    %Now realign and trim all data to the models
    fprintf('Aligning all data to the models\n');
    if ~isempty(fhM)
        Z2fhM = alignpulses2model(Z, fhM);
        Z2fhM = scaleZ2M(Z2fhM, fhM, scaleType);
        %trim data and model down to relevant parts (no padding) for first harmonic
        if max_length ~= lengthfhM
            startM = find(abs(fhM) > 0, 1,'first');
            finishM = find(abs(fhM) > 0, 1,'last');
            fhM = fhM(startM:finishM);
            S_fhM = S_fhM(:, startM:finishM);
            Z2fhM = Z2fhM(:, startM:finishM);
        end
    end
    clear startM finishM
    if ~isempty(shM)
        Z2shM = alignpulses2model(Z, shM);
        Z2shM = scaleZ2M(Z2shM, shM, scaleType);
        %trim data and model down to relevant parts (no padding) for second harmonic
        if max_length ~= lengthshM %max_length is from data (Z)
            startM = find(abs(shM) > 0, 1, 'first');
            finishM = find(abs(shM) > 0, 1, 'last');
            shM = shM(startM:finishM);
            S_shM = S_shM(:, startM:finishM);
            Z2shM = Z2shM(:, startM:finishM);
        end
    end
    clear startM finishM

    % calculate likelihood of data under each model
    if ~isempty(fhM)
        % Get standard deviation at each point
        S_ar_fh = repmat(S_fhM, size(Z2fhM, 1), 1);
        %calculate likelihood of data under each model
        fhM_ar = repmat(fhM, size(Z2fhM, 1), 1);
        LL_fhM = nansum(log10(normpdf(Z2fhM, fhM_ar, S_ar_fh)), 2);
        LL_0_fhpdf = nansum(log10(normpdf(Z2fhM, 0, S_ar_fh)), 2);
        LLR_fh = LL_fhM - LL_0_fhpdf;
    else
        LLR_fh = [];
    end

    if ~isempty(shM)
        % Get standard deviation at each point
        S_ar_sh = repmat(S_shM, size(Z2shM, 1), 1);
        % calculate likelihood of data under each model
        shM_ar = repmat(shM, size(Z2shM, 1), 1);
        LL_shM = nansum(log10(normpdf(Z2shM, shM_ar, S_ar_sh)), 2);
        LL_0_shpdf = nansum(log10(normpdf(Z2shM, 0, S_ar_sh)), 2);
        LLR_sh = LL_shM - LL_0_shpdf;
    else
        LLR_sh = [];
    end

    % pulse model + new-align-scaled-trimmed-data
    pulse_model.fhM = fhM;
    pulse_model.shM = shM;
    pulse_model.Z2fhM = Z2fhM; %aligned all pulses to first harmonic model
    pulse_model.Z2shM = Z2shM; %aligned all pulses to second harmonic model

    % generating pLik
    pLik.ModIdx = zeros(size(Z, 1), 1);
    if ~isempty(fhM) && ~isempty(shM)
        [best_LLR, mIdx] = max([LLR_fh, LLR_sh], [], 2);
        pLik.ModIdx(:, 1) = mIdx;
    elseif ~isempty(fhM)
       best_LLR = LLR_fh; pLik.ModIdx(:, 1) = 1;
    elseif ~isempty(shM)
       best_LLR = LLR_sh; pLik.ModIdx(:, 1) = 1;
    end

    pLik.LLR_best = best_LLR;
    pLik.LLR_fh = LLR_fh;
    pLik.LLR_sh = LLR_sh;
else
    pulse_model.fhM = [];
    pulse_model.shM = [];
    pulse_model.Z2fhM = []; %aligned all pulses to first harmonic model
    pulse_model.Z2shM = []; %aligned all pulses to second harmonic model
    pLik.LLR_best = zeros([], 1);
    pLik.LLR_fh = zeros([], 1);
    pLik.LLR_sh = zeros([], 1);
    pLik.ModIdx = zeros([], 2);
end

function [fZ, scale] = scaleZ2M(Z, M, type)
%fZ = scaleZ2M(Z,M,type)
%   type is an integer indicating what scalings should be applied.
%   Default value is 1, if type is not specified as an argument. 
%   Any value other than 1 indicates that only the first rescaling is 
%   used, and not the second.
if ~exist('type', 'var'); type = 1; end
[n_samples, total_length] = size(Z);
%rescale data to mean
%equivalent to following, on whole array
%a = mean(M.*Z(n,:))/mean(Z(n,:).^2);
%Z(n,:) = a*Z(n,:);
%% Scaling #1
Ma = repmat(M, n_samples, 1); % make as many copies of model as the # of signals.
num = mean(Ma'.*Z', 1); den = mean(Z'.^2, 1);
a = num./den;
ar = repmat(a', 1, total_length);
Z = ar.* Z;
scale = 1;
%% Scaling #2 - legacy code - do not use
if type == 1
    MZ = mean(Z,1); scale = max([(MZ/M) 0.30]); % 
end
fZ = Z/scale;

function fZ = alignpulses2model(Z, M)
% Horizontal alignment of the signal to the model.
[n_samples, total_length] = size(Z);
for n = 1:n_samples
    C = xcorr(M, Z(n,:),'unbiased');
    [~, tpeak] = max(C);
    tpeak = tpeak - total_length;
    Z(n,:) = circshift(Z(n,:),[1 tpeak]);
end
fZ = Z;