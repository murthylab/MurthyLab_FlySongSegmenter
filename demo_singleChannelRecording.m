cc()
addpath(genpath('src'))
%% load data
load('dat/PS_20130625111709_ch10.mat');
load('dat/PS_20130625111709_ch10manual.mat', 'pulseTimes');
minPulseTime = 30;%s
maxPulseTime = 120;%s

pulseTimes(pulseTimes<minPulseTime | pulseTimes>maxPulseTime)=[]; % cut to speedup demo
pulseTimesManual = pulseTimes-minPulseTime;
Fs = Data.fs; %Hz
%%
recording = Data.d(minPulseTime*Fs:maxPulseTime*Fs,:);
channels = size(recording, 2);
T = (1:size(recording,1))/Fs;
subplot(311)
plot(T, recording)
xlabel('time [ms]')
ylabel('voltage [V]')
axis('tight')
drawnow
%% process each channel
for chn = 1:channels
   [sInf(chn).nLevel, sInf(chn).winSine, sInf(chn).pulseInfo, sInf(chn).pulseInfo2, sInf(chn).pcndInfo] = ...
      segmentSong(recording(:,chn), 'params.m');
end
%% post process
bufferLen =  2e3;
noiseSample = findNoise(recording, bufferLen);
oneSong = recording; % same for single-channel data
[sInf, pInf, wInf, bInf, Song] = postProcessSegmentation(sInf, recording, oneSong, noiseSample);
pulseTimesAutomatic = pInf.wc/Fs;%s
%% plot results
clf
T = (1:length(oneSong))/Fs;
subplot(3,1,1:2)
plot(T, oneSong)
hold on
plot(pulseTimesAutomatic, ones(size(pulseTimesAutomatic))/5,'.','MarkerSize', 12)
plot(pulseTimesManual, ones(size(pulseTimesManual))/5+0.05, '.', 'MarkerSize',12)
axis('tight')
set(gca,'YLim', [-0.5 0.5], 'YTick', [0.2 0.25], 'YTickLabel', {'automatic', 'manual'})

subplot(313)
plot(T, bInf.Mask);
axis('tight')
set(gca, 'YTick', 0:2, 'YTickLabel', {'silence/noise', 'sine','pulse'})

linkaxes(gcas, 'x')
%% identify common events
% pulseTimes = [pulseTimesManual; pulseTimesAutomatic];
% pulseGroup = [ones(size(pulseTimesManual)); 2*ones(size(pInf.wc))];
% pulseGroupLabel = {'manual', 'automatic'}

% % pool all pulses
% sortIdx = argsort(pulseTimes);
% pulseTimes = pulseTimes(sortIdx);
% pulseGroup = pulseGroup(sortIdx);

% pulseTimes = pulseTimes(1:500);
% pulseGroup = pulseGroup(1:500);
% % go through pulses and pulses within +/- jitter ms
% jitter = 5/1000;%ms
% cnt = 1;
% pulseId = nan(size(pulseTimes));
% pulseId(1) = 1;
% for pul = 2:length(pulseTimes)
%    if pulseTimes(pul)-jitter>pulseTimes(pul-1)
%       cnt = cnt+1;
%    end
%    pulseId(pul) = cnt;
% end

% ids = max(pulseId);
% eventMat = zeros(ids, max(pulseGroup));
% for g = 1:max(pulseGroup)
%    eventMat(pulseId(pulseGroup==g), g) = 1;
% end
tolerance = 5/1000;%s
[confMat, eventMat] = idPulses(pulseTimesManual, pulseTimesAutomatic, tolerance)
%%
fprintf('\n')
fprintf('detected %d/%d pulses\n', length(pulseTimesAutomatic), length(pulseTimesManual))
fprintf('   true positives %d (p=%1.2f)\n',  sum(eventMat(:,1)==1 & eventMat(:,2)==1), mean(eventMat(:,1)==1 & eventMat(:,2)==1))
fprintf('   false negatives %d (p=%1.2f)\n', sum(eventMat(:,1)==1 & eventMat(:,2)==0), mean(eventMat(:,1)==1 & eventMat(:,2)==0))
fprintf('   true negative %d (p=%1.2f) (not really meaningful in this context)\n',   sum(eventMat(:,1)==0 & eventMat(:,2)==0), mean(eventMat(:,1)==0 & eventMat(:,2)==0))
fprintf('   false positives %d (p=%1.2f)\n', sum(eventMat(:,1)==0 & eventMat(:,2)==1), mean(eventMat(:,1)==0 & eventMat(:,2)==1))

