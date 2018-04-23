addpath(genpath('src'))
cc()
%% load data
load('dat/161118_1541bin');
load('dat/161118_1541bin_manual', 'pulseTimes');
pulseTimesManual = pulseTimes;
Fs = 10000;
%%
recording = double(data(6.4e5:9e5,1:9))./dataScalingFactor;
channels = size(recording, 2);
T = (1:size(recording,1))/Fs;
figure('Name', 'song segmentation')
clf
subplot(311)
plot(T, recording)
hold on
plot(pulseTimesManual, ones(size(pulseTimesManual))/2, '.', 'MarkerSize', 12)
xlabel('time [s]')
ylabel('voltage [V]')
axis('tight')
drawnow
%% segment song in each channel
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
pulseTimesAutomatic = pInf.wc/Fs;
%%
subplot(312)
plot(T, oneSong)
hold on
plot(pulseTimesAutomatic, ones(size(pulseTimesAutomatic))/5,'.','MarkerSize', 12)
plot(pulseTimesManual, ones(size(pulseTimesManual))/5+0.05, '.', 'MarkerSize',12)
set(gca,'YLim', [-0.4 0.4], 'YTick', [0.2 0.25], 'YTickLabel', {'automatic', 'manual'})

subplot(313)
plot(T, bInf.Mask);
set(gca, 'YTick', 0:2, 'YTickLabel', {'silence/noise', 'pulse','sine'})
axis(gcas, 'tight')
linkaxes(gcas, 'x')

%% compare automatic and manual segmentation
% identify pulse times in the manual and automatic data correspoding to the
% same pulse in the recording with a jitte of `tolerance` seconds
tolerance = 5/1000;%s
[confMat, eventMat] = idPulses(pulseTimesManual, pulseTimesAutomatic, tolerance);
confMatNorm = confMat./sum(confMat,1);

fprintf('\n')
fprintf('Detected %d/%d pulses:\n', sum(eventMat(:,2)==1), sum(eventMat(:,1)==1))
fprintf('   - true positives %d (p=%1.2f of all manually annotated pulses)\n',  sum(eventMat(:,1)==1 & eventMat(:,2)==1), sum(eventMat(:,1)==1 & eventMat(:,2)==1)./sum(eventMat(:,1)==1))
fprintf('   - false negatives %d (p=%1.2f of all manually annotated pulses)\n', sum(eventMat(:,1)==1 & eventMat(:,2)==0), sum(eventMat(:,1)==1 & eventMat(:,2)==0)./sum(eventMat(:,1)==1))
fprintf('   - false positives %d (p=%1.2f of all automatically called pulses)\n', sum(eventMat(:,1)==0 & eventMat(:,2)==1), sum(eventMat(:,1)==0 & eventMat(:,2)==1)./sum(eventMat(:,2)==1))
fprintf('   - true negatives %d (p=%1.2f, not meaningful in this context)\n',   sum(eventMat(:,1)==0 & eventMat(:,2)==0), sum(eventMat(:,1)==0 & eventMat(:,2)==0)./sum(eventMat(:,1)==0))


%% pulse type classification
pulsesNorm = normalizePulses(double(pInf.pSec)/1000);           % normalize pulses
pulseLabels = classifyPulses(pulsesNorm);       % classify pulses - 0=Pfast, 1=Pslow
fprintf('%d/%d Pfast, %d/%d Pslow.\n', sum(pulseLabels==0), length(pulseLabels), sum(pulseLabels==1), length(pulseLabels))

figure('Name', 'pulse type classification')
clf
subplot(311)
T = (1:size(pInf.pSec,2))/10;%ms
plot(T, double(pInf.pSec')/1000, 'Color', [0 0 0 0.2])
title('raw pulses')
subplot(312)
plot(T, pulsesNorm', 'Color', [0 0 0 0.2])
title('normalized pulses')
subplot(313)
hF = plot(T, pulsesNorm(pulseLabels==1,:)', 'Color', [1 0 0 0.2]);
hold on
hS = plot(T, pulsesNorm(pulseLabels==0,:)', 'Color', [0 0 0 0.2]);
hL = legend([hF(1) hS(1)], {' P_{fast} (red)', 'P_{slow} (black)'}, 'Box', 'off');
title('labeled pulses')
axis(gcas, 'tight')
set(gcas, 'Box', 'off', 'Color', 'none')

