addpath(genpath('src'))
cc()
%% load data
% CantonS recording from Stern (2014)
load('dat/PS_20130625111709_ch3.mat')
% hand-annotated pulse times for that recording from Kyriacou et al. (2017)
load('dat/PS_20130625111709_ch3manual.mat', 'pulseTimes')

% cut recording to the part that is hand-annotated and shorten to speed up
% processing for this demo
minPulseTime = 180;%s
maxPulseTime = 270;%max(pulseTimes);
pulseTimes(pulseTimes<minPulseTime | pulseTimes>maxPulseTime)=[];
pulseTimesManual = pulseTimes-minPulseTime;

Fs = Data.fs; %Hz
recording = Data.d(minPulseTime*Fs:maxPulseTime*Fs,:);
channels = size(recording, 2); % recording is time x channels

%% plot raw data
T = (1:size(recording,1))/Fs;
clf
subplot(3,1,1:2)
plot(T, recording)
hold on
plot(pulseTimesManual, ones(size(pulseTimesManual))/5+0.05, '.', 'MarkerSize',12)
xlabel('time [ms]')
ylabel('voltage [V]')
title('recording trace with manually annotated pulses')
axis('tight')
drawnow

%% process recording - detect sine and pulse
[sInf(chn).nLevel, sInf(chn).winSine, sInf(chn).pulseInfo, sInf(chn).pulseInfo2, sInf(chn).pcndInfo] = ...
   segmentSong(recording(:,chn), 'params.m');

%% post process
% automatically identify recording of duration `bufferLen` samples that does not
% contain song for estimating the noise floor in the recording
bufferLen =  2e3;
noiseSample = findNoise(recording, bufferLen);

% clean up pulses, identify bouts etc.
oneSong = recording; % same for single-channel data
[sInf, pInf, wInf, bInf, Song] = postProcessSegmentation(sInf, recording, oneSong, noiseSample);
pulseTimesAutomatic = pInf.wc/Fs;%s

%% plot results
clf
T = (1:length(oneSong))/Fs;
subplot(3,1,1:2)
plot(T, oneSong, 'k')
hold on
plot(pulseTimesAutomatic, ones(size(pulseTimesAutomatic))/5,'.','MarkerSize', 12)
plot(pulseTimesManual, ones(size(pulseTimesManual))/5+0.05, '.', 'MarkerSize',12)
axis('tight')
set(gca,'YLim', [-0.4 0.4], 'YTick', [0.2 0.25], 'YTickLabel', {'automatic', 'manual'})

subplot(313)
plot(T, bInf.Mask, 'k');
axis('tight')
set(gca, 'YTick', 0:2, 'YTickLabel', {'silence/noise', 'sine','pulse'})
xlabel('time [ms]')
set(gcas, 'Box','off','Color','none','TickDir','out')
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
