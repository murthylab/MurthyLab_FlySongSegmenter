function [confMat, eventMat, pulseId, pulseTimes, pulseGroup] = idPulses(pulseTimesA, pulseTimesB, tolerance)
% [eventMat, pulseId, pulseTimes, pulseGroup] = idPulses(pulseTimesA, pulseTimesB, tolerance)

if ~exist('tolerance','var') | isempty(tolerance)
	tolerance = 5/1000
end

% pool and sort all pulses
pulseTimes = [pulseTimesA; pulseTimesAutomatic];
pulseGroup = [ones(size(pulseTimesA)); 2*ones(size(pInf.wc))];

sortIdx = argsort(pulseTimes);
pulseTimes = pulseTimes(sortIdx);
pulseGroup = pulseGroup(sortIdx);

% go through pulses and pulses within +/- jitter ms
cnt = 1;
pulseId = nan(size(pulseTimes));
pulseId(1) = 1;
for pul = 2:length(pulseTimes)
   if pulseTimes(pul)-tolerance>pulseTimes(pul-1)
      cnt = cnt+1;
   end
   pulseId(pul) = cnt;
end

ids = max(pulseId);
eventMat = zeros(ids, max(pulseGroup));
for g = 1:max(pulseGroup)
   eventMat(pulseId(pulseGroup==g), g) = 1;
end

confMat = confusionMat(eventMat(:,1), eventMat(:,2))