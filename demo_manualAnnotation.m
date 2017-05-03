addpath(genpath('src'))
cc()
%% load data
load('dat/161118_1541bin');
Fs = 10000;
recording = double(data(6.4e5:9e5,1:9))./dataScalingFactor;
%%
% hitting the save button will save the manual annotation to 'workspace_byhand.mat'
% need to maximize window to reveal GUI elements (at least on OSX)
FlySongSegmenterByHand(recording, Fs)

%% plot results
load('workspace_byhand.mat')
delete('workspace_byhand.mat') % keep folder clean
pulseTimesManual = PULSE(:,2);
channels = size(recording, 2);
T = (1:size(recording,1))/Fs;
clf
plot(T, recording)
hold on
plot(pulseTimesManual, ones(size(pulseTimesManual))/2, '.', 'MarkerSize', 12)
xlabel('time [s]')
ylabel('voltage [V]')
axis('tight')
drawnow
