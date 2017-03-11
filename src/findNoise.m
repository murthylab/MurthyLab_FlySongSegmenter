function [noiseSample, noiseIdxRange] = findNoise(data, bufferLen)
% [mNoiseSamp, noiseIdxRange] = findNoise(data, bufferLen)
%
% uses a heuristic to detect a stretch of length bufferLen w/o song
% on all channels - may not work well for very dense recordings.
% ARGS:
%  data - time x channels
%  bufferLen - duration of the noise to find in samples (OPTIONAL - defaults to 10000);
% RETURNS:
%  noiseSample   - bufferLen x channels 
%  noiseIdxRange - start and end samples

if ~exist('bufferLen','var')
    bufferLen=10000; % samples
end
nChannels = size(data,2);

% compute smoothed env
env = abs(hilbert(data));
env = mapFun(@smooth,env,1000);

% compute power for different time windows 
for chn = 1:nChannels % per channel
   [envM, ~] = buffer(env(bufferLen:end-bufferLen,chn), bufferLen); 
   envVar(:,chn) = var(envM);
end

% select time window that has no song in it
% this is based on heuristics and only works for sufficiently sparse recordings
envVarMax = max(log2(envVar),[],2);% max var across channels 
minIdx = argmin(envVarMax);   % take time in wich the max is minimal

% outputs
% noiseIdxRange = bufferLen + (minIdx-1)*bufferLen + [1 bufferLen];
noiseIdxRange = bufferLen + (minIdx-1)*bufferLen + [1 bufferLen];
noiseSample = data(noiseIdxRange(1):noiseIdxRange(2),:);
