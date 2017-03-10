cc()
addpath(genpath('src'))
%% load data
load('dat/161118_1541bin');
Fs = 10000;
%%
recording = double(data(6.4e5:9e5,1:9))./dataScalingFactor;
channels = size(recording, 2);
T = (1:size(recording,1))/Fs;
subplot(311)
plot(T, recording)
xlabel('time [ms]')
ylabel('voltage [V]')
axis('tight')
drawnow
%% preprocess (if necessary)
for chn = 1:channels
   fprintf('segmenting channel %d.\n', chn)
   [sInf(chn).nLevel, sInf(chn).winSine, sInf(chn).pulseInfo, sInf(chn).pulseInfo2, sInf(chn).pcndInfo] = ...
      segmentSong(recording(:,chn), 'params.m');
end
%% post process
bufferLen =  2e3; % samples
noiseSample = findNoise(recording, bufferLen);
oneSong = mergeChannels(recording)';
[sInf, pInf, wInf, bInf, Song] = postProcessSegmentation(sInf, recording, oneSong, noiseSample);
%%
subplot(312)
plot(T, oneSong)
hold on
plot(pInf.wc/Fs, ones(size(pInf.wc))/5,'.','MarkerSize', 12)
subplot(313)
plot(T, bInf.Mask);
set(gca, 'YTick', 0:2, 'YTickLabel', {'silence/noise', 'pulse','sine'})
axis(gcas, 'tight')
linkaxes(gcas, 'x')

