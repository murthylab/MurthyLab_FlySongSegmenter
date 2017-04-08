% no - do this based on parameter
% assemble everything into struct
cDir = pwd;

% decides which model to used based on the folder name
if ~isempty(strfind(cDir,'VT40556')) % params for pIP10   
    % ******************* Pulse related variables *******************
    if ~isempty(strfind(cDir,'LG'))
        PulseModel2Use = 'VT40556-TRPA1LG'; % pulse model used for pLik
    elseif ~isempty(strfind(cDir,'Free'))
        PulseModel2Use = 'VT40556-TRPA1Free'; % pulse model used for pLik
    else
        PulseModel2Use = 'VT40556-Reach'; % pulse model used for pLik
    end
    minPulseNum = 2; % min number of pulses per recording;
    pLik_Thrs = 0; % min pLik, likelihood to model pulse, to accept "pLikThresholding"
    pLik_Thrs_sPul = 16; % min pLik for single pulse
    pLik_Thrs_nPul = [-120 -20]; % range of pLik to grow bouts
    pLik_scaletype = 0; % which kind of scaling you do to traces for pLik calcaulation
    ampR_Thrs = 1.744; % min of relative amplitude
    ampR_Thrs_sPul = 2.3; % min ampR for single pulse

    % ******************* sine related variables used by "combineSineSong" *******************
    sine_nThrs = 0.15; % relative to the number of timpoints that has sine > ths
    sine_LThrs = 700; % min length of sine epocs (default 70 ms)

    % ******************* bout related variables *******************
    % Variables used to create bouts within "findSongBouts"
    TimeThrs = 3000; % delete bouts / pulses close to start/end of recording (default: 35000)
    buff = [TimeThrs TimeThrs]; % buffer to add to each bout found
    s2s_disThrs = 2000; % song 2 song distance threshold for bout definition
    p2p_disThrs = 10^4; % max pulse to pulse / bout distance before deleting it
    p2p_growbout_disThrs = 1100; % pulse 2 pulse distance for growing bouts
    % variables used to delete bouts with low pulse/sine
    minPn = [1 1 1 0 1 1]; % unit: pulse number % [2 2 2 0 2 2]
    minSl = [1 2 1 3 5]; % unit: hundreds of ms

    % ******************* secondary parameters to look at *******************
    hpWi_Thrs = [84 250 480]; % threshold of hpWi
    min_max_fcmx = [125 225]; % wavelet frequency at peak 
    
    % ******************* Implement protocol specific switch *******************
    noiseThsGate = 1;
    if ~isempty(strfind(cDir,'LG')) || ~isempty(strfind(cDir,'Free'))
        fprintf('Overlooking noise level (rej vs true pulses during stim vs non-stim)\n')
        noiseThsGate = 0;
    end
    
elseif ~isempty(strfind(cDir,'NP2631')) && isempty(strfind(cDir,'NP2631-FruFLP'))
    % ******************* Pulse related variables *******************
    if ~isempty(strfind(cDir,'LG'))
        PulseModel2Use = 'NP2631-TRPA1LG'; % pulse model used for pLik
    elseif ~isempty(strfind(cDir,'Free'))
        PulseModel2Use = 'NP2631-TRPA1Free'; % pulse model used for pLik
    else
        PulseModel2Use = 'NM91'; %'VT40556-Reach'; % pulse model used for pLik
    end
    minPulseNum = 2; % min number of pulses per recording;
    pLik_Thrs = 0; % min pLik, likelihood to model pulse, to accept "pLikThresholding"
    pLik_Thrs_sPul = 16; % min pLik for single pulse
    pLik_Thrs_nPul = [-120 -20]; % range of pLik to grow bouts
    pLik_scaletype = 0; % which kind of scaling you do to traces for pLik calcaulation
    ampR_Thrs = 1.744; % min of relative amplitude
    ampR_Thrs_sPul = 2.3; % min ampR for single pulse

    % ******************* sine related variables used by "combineSineSong" *******************
    sine_nThrs = 0.15; % relative to the number of timpoints that has sine > ths
    sine_LThrs = 700; % min length of sine epocs (default 70 ms)

    % ******************* bout related variables *******************
    % Variables used to create bouts within "findSongBouts"
    TimeThrs = 3000; % delete bouts / pulses close to start/end of recording (default: 35000)
    buff = [TimeThrs TimeThrs]; % buffer to add to each bout found
    s2s_disThrs = 2000; % song 2 song distance threshold for bout definition
    p2p_disThrs = 10^4; % max pulse to pulse / bout distance before deleting it
    p2p_growbout_disThrs = 1100; % pulse 2 pulse distance for growing bouts
    % variables used to delete bouts with low pulse/sine
    minPn = [1 1 1 0 1 1]; % unit: pulse number % [2 2 2 0 2 2]
    minSl = [1 2 1 3 5]; % unit: hundreds of ms

    % ******************* secondary parameters to look at *******************
    hpWi_Thrs = [84 250 480]; % threshold of hpWi
    min_max_fcmx = [125 225]; % wavelet frequency at peak 
    
    % ******************* Implement protocol specific switch *******************
    noiseThsGate = 0;
    
elseif ~isempty(strfind(cDir,'VT57239-FruFLP'))
    % ******************* Pulse related variables *******************
    PulseModel2Use = {'NM91', 'VT57239-FruFLP_m1', 'VT57239-FruFLP_m2', 'VT57239-FruFLP_m3'}; %'VT40556-Reach' 
    minPulseNum = 2; % min number of pulses per recording;
    pLik_Thrs = 0; % min pLik, likelihood to model pulse, to accept "pLikThresholding"
    pLik_Thrs_sPul = 16; % min pLik for single pulse
    pLik_Thrs_nPul = [-120 -20]; % range of pLik to grow bouts
    pLik_scaletype = 0; % which kind of scaling you do to traces for pLik calculation
    ampR_Thrs = 2.3; % min of relative amplitude
    ampR_Thrs_sPul = 3; % min ampR for single pulse

    % ******************* sine related variables used by "combineSineSong" *******************
    sine_nThrs = 0.15; % relative to the number of timpoints that has sine > ths
    sine_LThrs = 700; % min length of sine epocs (default 70 ms)

    % ******************* bout related variables *******************
    % Variables used to create bouts within "findSongBouts"
    TimeThrs = 3000; % delete bouts / pulses close to start/end of recording (default: 35000)
    buff = [TimeThrs TimeThrs]; % buffer to add to each bout found
    s2s_disThrs = 2000; % song 2 song distance threshold for bout definition
    p2p_disThrs = 10^4; % max pulse to pulse / bout distance before deleting it
    p2p_growbout_disThrs = 1100; % pulse 2 pulse distance for growing bouts
    % variables used to delete bouts with low pulse/sine
    minPn = [1 1 1 0 1 1]; % unit: pulse number % [2 2 2 0 2 2]
    minSl = [1 2 1 3 5]; % unit: hundreds of ms

    % ******************* secondary parameters to look at *******************
    hpWi_Thrs = [84 250 480]; % threshold of hpWi
    min_max_fcmx = [125 225]; % wavelet frequency at peak 
    
    % ******************* Implement protocol specific switch *******************
    noiseThsGate = 0;
    
elseif ~isempty(strfind(cDir,'NP2631-FruFLP')) || ~isempty(strfind(cDir,'P1-Split'))% params for P1
    % ******************* Pulse related variables *******************
    PulseModel2Use = 'VT40556-Reach'; % pulse model used for pLik
    minPulseNum = 2; % min number of pulses per recording;
    pLik_Thrs = 0; % min pLik, likelihood to model pulse, to accept "pLikThresholding"
    pLik_Thrs_sPul = 16; % min pLik for single pulse
    pLik_Thrs_nPul = [-120 -20]; % range of pLik to grow bouts
    pLik_scaletype = 0; % which kind of scaling you do to traces for pLik calcaulation
    ampR_Thrs = 1.744; % min of relative amplitude
    ampR_Thrs_sPul = 2.3; % min ampR for single pulse

    % ******************* sine related variables used by "combineSineSong" *******************
    sine_nThrs = 0.15; % relative to the number of timpoints that has sine > ths
    sine_LThrs = 700; % min length of sine epocs (default 70 ms)

    % ******************* bout related variables *******************
    % Variables used to create bouts within "findSongBouts"
    TimeThrs = 3000; % delete bouts / pulses close to start/end of recording (default: 35000)
    buff = [TimeThrs TimeThrs]; % buffer to add to each bout found
    s2s_disThrs = 2000; % song 2 song distance threshold for bout definition
    p2p_disThrs = 10^4; % max pulse to pulse / bout distance before deleting it
    p2p_growbout_disThrs = 1100; % pulse 2 pulse distance for growing bouts
    % variables used to delete bouts with low pulse/sine
    minPn = [1 1 1 0 1 1]; % unit: pulse number % [2 2 2 0 2 2]
    minSl = [1 2 1 3 5]; % unit: hundreds of ms

    % ******************* secondary parameters to look at *******************
    hpWi_Thrs = [84 250 480]; % threshold of hpWi
    min_max_fcmx = [125 225]; % wavelet frequency at peak 
    
    % ******************* Implement protocol specific switch *******************
    noiseThsGate = 0;
    if ~isempty(strfind(cDir,'P1-Split'))
        noiseThsGate = 0;
    end
end