cc()
addpath(genpath('src'))
%% load data
load('dat/PS_20130625111709_ch10.mat');
Fs = Data.fs; %Hz
%%
recording = Data.d(1:30e4,:);
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
%% plot results
T = (1:length(oneSong))/Fs;
subplot(3,1,1:2)
plot(T, oneSong)
hold on
plot(pInf.wc/Fs, ones(size(pInf.wc))/5,'.','MarkerSize', 12)
subplot(313)
plot(T, bInf.Mask);
set(gca, 'YTick', 0:2, 'YTickLabel', {'silence/noise', 'sine','pulse'})
axis(gcas,'tight')
linkaxes(gcas, 'x')


