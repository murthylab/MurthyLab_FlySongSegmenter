function pulses = normalizePulses(pulses, signWindow, smoothWindow, upsamplingFactor, restoreFlag, centerFlag, shiftWindow)
% pulses = normalizePulses(pulses, signWindow=[10 0], smoothWindow=15, restoreFlag=true)
% normalizes pulses by 1. dividing by norm, 2. aliging to peak energy (envelope), and 3. flipping sign
%
% ARGS
%  pulses - Npulses x Nsamples matrix of pulse shapes
%  signWindow   - n samples (preceding peak) to use for sign flipping
%  smoothWindow - n samples with which to smooth pulse for RMS envelope estimation
%  upsamplingFactor - for upsampling pulses - all other params (signWindow, smoothWindow) 
%                     refer to the upsampled scale
%  restoreFlag - downsample upsampled and aligned pulses or return upsampled pulses
%  centerFlag  - 
% RETURNS
%  normalized pulses

if ~exist('signWindow','var') || isempty(signWindow)
   signWindow = [10 0];
end
if ~exist('smoothWindow', 'var') || isempty(smoothWindow)
   smoothWindow = 15;
end
if ~exist('upsamplingFactor', 'var') || isempty(upsamplingFactor)
   upsamplingFactor = 1;
end
if ~exist('restoreFlag', 'var') || isempty(restoreFlag)
   restoreFlag = 1;
end
if ~exist('centerFlag', 'var') || isempty(centerFlag)
   centerFlag = 1;
end
if ~exist('shiftWindow','var') || isempty(shiftWindow)
   shiftWindow = [];
end


pulses = double(pulses);

% upsample if requested
if upsamplingFactor>1
   try
      pulses = resample(pulses', upsamplingFactor, 1, 100)';
   catch ME
      disp(ME.getReport())
   end
end
[nPulses, pulseDur] = size(pulses);
if isempty(shiftWindow)
   shiftIdx = 1:pulseDur;
else
   shiftIdx = ceil(pulseDur/2) + [-shiftWindow:shiftWindow];
end
%% Center on max power and flip if negative lobe to left of max power
for ii = 1:nPulses
   oldPulse = pulses(ii,:);
   oldPulse = oldPulse./norm(oldPulse);                     % scale to unit norm
   
   % center of mass
   switch centerFlag
      case 1
         mIdx = argmax(smooth(oldPulse(shiftIdx).^2, smoothWindow));       % get max power then center on this
      case 2
         mIdx = argmax(smooth(abs(hilbert(oldPulse(shiftIdx))), smoothWindow));       % get max power then center on this
      case 3
         mIdx = argmax(smooth(abs(oldPulse(shiftIdx)), smoothWindow));       % get max power then center on this
      case 4
         mIdx = round(centerOfMass((1:pulseDur), smooth(oldPulse(shiftIdx).^2, smoothWindow)'));       % get max power then center on this
      case 5
         mIdx = round(centerOfMass((1:pulseDur), smooth(abs(hilbert(oldPulse(shiftIdx))), smoothWindow)'));       % get max power then center on this
      case 6
         mIdx = round(centerOfMass((1:pulseDur), smooth(abs(oldPulse(shiftIdx)), smoothWindow)'));       % get max power then center on this
      otherwise
         error()
   end
   % max energy
   lPad = pulseDur - mIdx - min(shiftIdx);  %i.e. ((total_length/2) - C)    
   rPad = 2*pulseDur - pulseDur - lPad;
   newPulse = [zeros(1,lPad) oldPulse zeros(1, rPad)];
   if mean(newPulse(pulseDur-signWindow(1):pulseDur-signWindow(2)))<0
      newPulse = -newPulse;                                 % flip pulse if it starts neg.
   end
   pulses(ii,:) = newPulse(floor(pulseDur/2)+(1:pulseDur)); % cut normalized pulse to original duration
end

%% restore original pulse sampling
if upsamplingFactor>1 && restoreFlag==1
   try
      pulses = resample(pulses', 1, upsamplingFactor, 100)';
   catch ME
      disp(ME.getReport())
   end
end
