function [confMat, eventMat, pulseId, pulseTimes, pulseGroup] = idPulses(pulseTimesA, pulseTimesB, tolerance)
% [eventMat, pulseId, pulseTimes, pulseGroup] = idPulses(pulseTimesA, pulseTimesB, tolerance)

if ~exist('tolerance','var') || isempty(tolerance)
	tolerance = 5/1000;% 5ms
end

% pool and sort all pulses
pulseTimes = [pulseTimesA; pulseTimesB];
pulseGroup = [ones(size(pulseTimesA)); 2*ones(size(pulseTimesB))];

sortIdx = argsort(pulseTimes);
pulseTimes = pulseTimes(sortIdx);
pulseGroup = pulseGroup(sortIdx);

% go through pulses and assign pulses within tolerance to same event
cnt = 1;
pulseId = nan(size(pulseTimes));
pulseId(1) = 1;
for pul = 2:length(pulseTimes)
   if pulseTimes(pul)-tolerance>pulseTimes(pul-1)
      cnt = cnt+1;
   end
   pulseId(pul) = cnt;
end

% build binary event matrix for each set of pulse times
ids = max(pulseId);
eventMat = zeros(ids, max(pulseGroup));
for g = 1:max(pulseGroup)
   eventMat(pulseId(pulseGroup==g), g) = 1;
end

confMat = confusionmat(eventMat(:,1), eventMat(:,2));