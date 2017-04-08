function p = song_postProcessSegmentationParameters()

p.PulseModel2Use = 'NM91'; %'VT40556-Reach'; % pulse model used for pLik
p.minPulseNum = 2; % min number of pulses per recording;
p.pLik_Thrs = 0; % min pLik, likelihood to model pulse, to accept "pLikThresholding"
p.pLik_Thrs_sPul = 16; % min pLik for single pulse
p.pLik_Thrs_nPul = [-120 -20]; % range of pLik to grow bouts
p.pLik_scaletype = 0; % which kind of scaling you do to traces for pLik calcaulation
p.ampR_Thrs = 1.744; % min of relative amplitude
p.ampR_Thrs_sPul = 2.3; % min ampR for single pulse

% ******************* sine related variables used by "combineSineSong" *******************
p.sine_nThrs = 0.15; % relative to the number of timpoints that has sine > ths
p.sine_LThrs = 700; % min length of sine epocs (default 70 ms)

% ******************* bout related variables *******************
% Variables used to create bouts within "findSongBouts"
p.TimeThrs = 3000; % delete bouts / pulses close to start/end of recording (default: 35000)
p.buff = [p.TimeThrs p.TimeThrs]; % buffer to add to each bout found
p.s2s_disThrs = 2000; % song 2 song distance threshold for bout definition
p.p2p_disThrs = 10^4; % max pulse to pulse / bout distance before deleting it
p.p2p_growbout_disThrs = 1100; % pulse 2 pulse distance for growing bouts
% variables used to delete bouts with low pulse/sine
p.minPn = [1 1 1 0 1 1]; % unit: pulse number % [2 2 2 0 2 2]
p.minSl = [1 2 1 3 5]; % unit: hundreds of ms

% ******************* secondary parameters to look at *******************
p.hpWi_Thrs = [84 250 480]; % threshold of hpWi
p.min_max_fcmx = [125 225]; % wavelet frequency at peak
