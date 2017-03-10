function [pulse_model,Lik_pulse] = fitPulse2GivenModel(pulses, pMod)
% [pulse_model,Lik_pulse] = Z_2_pulse_model(pulse_model,new_pulses,[sample_freqs])
% USAGE
%  If the model and data were sampled at different sampling frequnecies, then enter this in sample_freqs
%  e.g. [1e4 4e4] as the freqs for the model and new_pulses
%  provide sample of pulses
% RETURNS
%  pulse model & std etc and Lik of individual pulses given the model

fhM = pMod.fhM;
shM = pMod.shM;
lengthfhM = length(fhM);
lengthshM = length(shM);
S_fhM = double(pMod.S_fhM);
S_shM = double(pMod.S_shM);
d = pulses;
if ~iscell(d);
   if size(d,1) == 251
      d = d';
   end
   d = mat2cell(d', 251, ones(size(d,1),1));
end

% grab samples, center, and pad
n_samples = length(d);
max_length = max(cellfun(@length,d));
total_length = 2* max_length;
Z = zeros(n_samples,total_length );
for n=1:n_samples;
   X = d{n};
   T = length(X);
   [~,C] = max(X);%get position of max power
   
   %center on max power
   left_pad = max_length - C; %i.e. ((total_length/2) - C)
   right_pad = total_length - T - left_pad;
   Z(n,:) = [zeros(left_pad,1); X ;zeros((right_pad),1)];
end

% pad models
fhPad = round(lengthfhM /2);
shPad = round(lengthshM /2);
fhM = [zeros(fhPad,1)', fhM , zeros(fhPad,1)'];
shM = [zeros(shPad,1)', shM , zeros(shPad,1)'];
S_fhM = [zeros(size(S_fhM,1),fhPad,1), S_fhM , zeros(size(S_fhM,1),fhPad,1)];
S_shM = [zeros(size(S_shM,1),shPad,1), S_shM , zeros(size(S_shM,1),shPad,1)];

%% double check to ensure lengths of model and data are equal
% first harmonic first
if ~isequal(total_length,length(fhM))
   diff = total_length - length(fhM);
   if diff > 1 %if data longer than model, add zeros(diff,1) to end of model
      fhM = [fhM, zeros(diff,1)'];
   else %add zeros to end of data
      Z = [Z,zeros(size(Z,1),-diff)];%take neg diff, because is neg
   end
end

% then second harmonic
if ~isequal(size(Z,2),length(shM))
   diff = size(Z,2) - length(shM);
   if diff > 1 % if data longer than model
      shM = [shM, zeros(diff,1)'];
   else % this is a very unlikely scenario, go ahead and trim second harmonic model is somehow longer than data
      if mod(diff,2) % if diff is odd
         left_trim =round(diff/2);
         right_trim = left_trim - 1 ;
         shM = shM(left_trim:end-right_trim);
      else
         trim = diff/2;
         shM= shM(trim:end-trim);
      end
   end
end


% align new data to old models
Z2fhM = alignpulses2model(Z,fhM);
Z2fhM = scaleZ2M(Z2fhM,fhM);
Z2shM = alignpulses2model(Z,shM);
Z2shM = scaleZ2M(Z2shM,shM);

% Generate phase reversed model
fhRM = -fhM;
Z2fhRM = alignpulses2model(Z,fhRM);
Z2fhRM = scaleZ2M(Z2fhRM,fhRM);

% Generate phase reversed second harmonic model
shRM = -shM;
Z2shRM = alignpulses2model(Z,shRM);
Z2shRM = scaleZ2M(Z2shRM,shRM);

chisq = nan(n_samples, 4);
for n=1:n_samples
   %calc chi-square for first model
   chisq(n,1) = ...
      mean((Z2fhM(n,:) - fhM).^2./var(Z2fhM(n,:)));
   % calc chi-square for second harmonic model
   chisq(n,2) = ...
      mean((Z2shM(n,:) - shM).^2./var(Z2shM(n,:)));
   % calc chi-square for reversed model
   chisq(n,3) = ...
      mean((Z2fhRM(n,:) - fhRM).^2./var(Z2fhRM(n,:)));
   % calc chi-square for reversed second harmonic model
   chisq(n,4) = ...
      mean((Z2shRM(n,:) - shRM).^2./var(Z2shRM(n,:)));
end
[~,best_chisqr_idx] = nanmin(chisq,[],2);

% flip data that fits a reversed model better (columns 3 or 4)
for n=1:n_samples
   if best_chisqr_idx(n) > 2
      Z(n,:) = -Z(n,:);
   end
end
% Now realign all data to the models
fprintf('Aligning all data to the models.\n');
Z2fhM = alignpulses2model(Z,fhM);
[Z2fhM, scaleFHM] = scaleZ2M(Z2fhM,fhM);
Z2shM = alignpulses2model(Z,shM);
[Z2shM, scaleSHM] = scaleZ2M(Z2shM,shM);

%trim data and model down to relevant parts (no padding) for first harmonic
%If length of original data is longer than model
%then trim data to length of original models
%if data longer than model, just trim to model
if max_length > lengthfhM %max_length is from data (Z)
   startM = find(abs(fhM>0),1,'first');
   finishM = find(abs(fhM>0),1,'last');
   fhM = fhM(startM:finishM);
   S_fhM = S_fhM(:,startM:finishM);
   
   Z2fhM = Z2fhM(:,startM:finishM);
   %if model longer than data, trim as follows
   %compare SE of Z at each point (from front and back) with deviation of fh model
   %start and stop when deviation exceeds SE of data
elseif max_length < lengthfhM
   %if model is longer than data, then trim model
   %compare SE at each point (from front and back) with deviation of fh model
   %start and stop when deviation exceeds SE of data
   S_Z = std(Z2fhM(Z2fhM ~= 0));%take only data that are not 0 (i.e. padding)
   SE_Z = S_Z/sqrt(n_samples);
   startZ = find((abs(mean(Z2fhM))>SE_Z),1,'first');
   finishZ = find((abs(mean(Z2fhM))>SE_Z),1,'last');
   fhM = fhM(startZ:finishZ);
   S_fhM = S_fhM(:,startZ:finishZ);
   Z2fhM = Z2fhM(:,startZ:finishZ);
end

%trim data and model down to relevant parts (no padding) for second harmonic
if max_length > lengthshM %max_length is from data (Z)
   startM = find(abs(shM>0),1,'first');
   finishM = find(abs(shM>0),1,'last');
   shM = shM(startM:finishM);
   S_shM = S_shM(:,startM:finishM);
   Z2shM = Z2shM(:,startM:finishM);
   %if model longer than data, trim as follows
   %compare SE of Z at each point (from front and back) with deviation of fh model
   %start and stop when deviation exceeds SE of data
elseif max_length < lengthshM
   %if model is longer than data, then trim model
   %compare SE at each point (from front and back) with deviation of fh model
   %start and stop when deviation exceeds SE of data
   S_Z = std(Z2shM(Z2shM ~= 0));%take only data that are not 0 (i.e. padding)
   SE_Z = S_Z/sqrt(n_samples);
   startZ = find((abs(mean(Z2shM))>SE_Z),1,'first');
   finishZ = find((abs(mean(Z2shM))>SE_Z),1,'last');
   shM = shM(startZ:finishZ);
   S_shM = S_shM(:,startZ:finishZ);
   Z2shM = Z2shM(:,startZ:finishZ);
end

% Get standard deviation at each point
% S_shM = std(shZ);

S_ar_fh = repmat(S_fhM,size(Z2fhM,1),1);
S_ar_sh = repmat(S_shM,size(Z2shM,1),1);

%% calculate likelihood of data under each model
fhM_ar = repmat(fhM,size(Z2fhM,1),1);
LL_fhM = nansum(log10(normpdf(Z2fhM,fhM_ar,S_ar_fh)),2);
LL_0_fhpdf = nansum(log10(normpdf(Z2fhM,0,S_ar_fh)),2);

shM_ar = repmat(shM,size(Z2shM,1),1);
LL_shM = nansum(log10(normpdf(Z2shM,shM_ar,S_ar_sh)),2);
LL_0_shpdf = nansum(log10(normpdf(Z2shM,0,S_ar_sh)),2);

LLR_fh = LL_fhM - LL_0_fhpdf;
LLR_sh = LL_shM - LL_0_shpdf;

% Take best LLR
best_LLR = max(LLR_fh,LLR_sh);
pulse_model.fhM = fhM;
pulse_model.shM = shM;
pulse_model.Z2fhM = Z2fhM;%aligned all pulses to first harmonic model
pulse_model.Z2shM = Z2shM;%aligned all pulses to first harmonic model
Lik_pulse.LLR_best = best_LLR;
Lik_pulse.LLR_fh = LLR_fh;
Lik_pulse.LLR_sh = LLR_sh;


