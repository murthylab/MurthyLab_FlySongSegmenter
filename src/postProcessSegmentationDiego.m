function	[sInf, pInf, wInf, bInf, Song, dInf] = postProcessSegmentationDiego(sInf, data, oneSong, mNoiseSamp, dataScalingFactor, lastSampleToProcess, pulseModelName, pulseModels, parameterFile)
%[sInf, pInf, wInf, bInf, Song] = postProcessSegmentation(sInf, data, oneSong, mNoiseSamp, dataScalingFactor, lastSampleToProcess, pulseModelName)

% OLD SIGNATURE: function [sInf, pInf, wInf, bInf, Song, dInf] = postProcessSegmentationDiego(sInf, oneSong, mNoiseSamp, data)
% general workflow
% get main pulse features to use:
% pulse centers (wc), mic detected (detM), center2surround ratio (c2SR)
% pulse trace per mic (pSec), pulse amplitude (pInf.mAmp.Int, pInf.mAmp.Mx, pInf.ampR)
% sine intervals (stEn), temporaly holds sine that pass amplitude threshold (aSin), mic detected (bSin)

% TODO:
% generate mNoiseSamp and oneSong on-the-fly from data if not supplied as arg


if ~exist('dataScalingFactor','var') || isempty(dataScalingFactor)
   dataScalingFactor = 1;
end
if ~exist('oneSong','var') || isempty(oneSong)
   warning('no merged channels provided - merging channels.')
   oneSong = mergeChannels(double(data)/dataScalingFactor);
end
if ~exist('mNoiseSamp','var') || isempty(mNoiseSamp)
   warning('no noise sample provided - estimating automatically.')
   mNoiseSamp = findNoise(double(data)/dataScalingFactor);
end
if ~exist('lastSampleToProcess','var') || isempty(lastSampleToProcess)
   lastSampleToProcess = size(data, 1); % default to last sample in recording
else
   warning('`lastSampleToProcess` not implemented yet - will process full song')
end
if ~exist('pulseModelName','var') || isempty(pulseModelName)
   pulseModelName = ''; % default to empty - will use NM91 model
end
if ~exist('parameterFile','var') || isempty(parameterFile)
   parameterFile = 'postProcessSegmentationDiego_Parameters'; % default to empty - will use NM91 model
end


% ******************************************************************************
% *********** load segmentation parameters from file ***************************
% ******************************************************************************
% TODO: specify name as parameter...
global p
p = eval(parameterFile);
% ******************************************************************************
% *********** Collecting main pulse and sine parameters (pInf, wInf) ***********
% ******************************************************************************
% 

% transpose song data (for compatibility with original postProcessing code)
data = double(data')/dataScalingFactor;
oneSong = oneSong';

% % getting pulse model
% load pulseModels;
% % get Idx for pulse models to use
% pM_Idx = [];
% if ~iscell(p.PulseModel2Use); p.PulseModel2Use = {p.PulseModel2Use}; end
% for pM_i = 1:numel(p.PulseModel2Use)
%    pM_Idx = [pM_Idx, find(cellfun(@(x) strcmp(p.PulseModel2Use{pM_i}, x), {pMod.fStr}') ~= 0)];
% end
% pMod = pMod(pM_Idx);

% load pulse models
if ~exist('pulseModels','var') || isempty(pulseModels)
   load('pulseModels', 'pMod');
else
   pMod = pulseModels;
end

if strfind(pulseModelName, 'NM91')
   pModIdx = cellfun(@(x) strcmpi('NM91', x), {pMod.fStr}');
else
   pModIdx = cellfun(@(x) strcmpi(pulseModelName, x), {pMod.fStr}');
end
if any(pModIdx)
   pMod = pMod(pModIdx);
else
   warning('did not find pulse model matching "%s". defaulting to NM91.', pulseModelName)
   pMod = pMod(cellfun(@(x) strcmpi('NM91', x), {pMod.fStr}'));
end


% getting noise amplitude
noiseOffset = mean(mNoiseSamp);
nAmp = abs(bsxfun(@minus, double(mNoiseSamp), noiseOffset))/1000;
for i = 1:size(nAmp, 2)
   nAmp(:, i) = runningExtreme(nAmp(:, i), 101, 'max');
end
nAmp = max(mean(nAmp));
wInf.nAmp = nAmp; clear nAmp

%% Combining mics: generating pInf, wInf var and proceed to pulse / sine removal
% getting pulse centers (pInf.Song, pInf.detM, pInf.wc, pInf.c2SR)
% fuse centers with distances lower than 5 ms, and distances
% within 5:12 ms are evaluated to find the center (average position).
if max(oneSong(:))>3; oneSong = double(oneSong)/1000; end
pInf.Song = oneSong;
[pInf.detM, pInf.wc, pInf.c2SR] = findUniPul(sInf, pInf.Song);

% delete pulses close to trial edges using "p.TimeThrs"
pInf = rmpulseclose2edge(pInf);

% updates wc based on frequency and also generates (pInf.dog, pInf.fcmx)
%pInf = fcmx_calculator(pInf, sInf, wInf, mNoiseSamp, noiseOffset, fName);

% center pulse on max power, center pSec and corrects wc
[pInf] = getpsec(pInf, sInf, wInf, noiseOffset, mNoiseSamp, data);
[pInf] = centerPulse(pInf);

% getting sine params (wInf.stEn, wInf.aSin, wInf.bSin)
% delete epocs of sine shorter than "p.sine_LThrs" (default 70 ms)
wInf = wInfGen(sInf, pInf, wInf);

% getting pSec (pulse trace per mic, wc+-25ms), their c2SR (pInf.c2SR, pInf.pSec) and amplitude per mic
% (pInf.mAmp.Int, integral, pInf.mAmp.Mx, max(abs(pulse)), pInf.ampR, relative amplitude to noise levels)
[pInf, wInf] = getp_Amp_c2SR(pInf, sInf, wInf ,noiseOffset, mNoiseSamp, data);
clear oneSong i mNoiseSamp data;

% getting hpWi per pulse (pInf.hpWi)
pInf = getp_hpWi_pulse2song(pInf);

% getting pulse2sine|pulse distance (pInf.proxSng) using p.p2p_disThrs
pInf = findSongNearPulse(pInf, wInf);

% generate pInf.exSong, used by remSinConPul
pInf.exSong = runningExtreme(abs(pInf.Song), 65, 'max');

% getting pInf.pLik, using (sInfFitPulse2Mod)
pInf.pLik = zeros([],1); pInf.Model = zeros([],1);
if numel(pInf.wc)> 0
   [pInf.Model, pInf.pLik] = sInfFitMultiPulseMod(pInf.pSec, pMod, p.pLik_scaletype);
end

fprintf(['Initial number of : pulses ',num2str(length(pInf.wc)),...
   ' , sine epocs ',num2str(length(wInf.stEn)/2),'\n'])

% generate pInf.boutTag, pInf.ppBout, pInf.boutLikMax, pInf.boutLikMin
[~, ~, pInf] = findSongBouts(pInf, wInf);

% local storage of pulse shapes for further offline inspection
dInf = pInf;

% *************************************************************
% ********************** Removing Pulses **********************
% *************************************************************

% remove pulses using a combination of
% 1) pLik to model (pLik >= 0)
% 2) distance to neighbouring pulse /sine, proxSng > 10^4
% 3) amplitude ampR (ampR > 1.744)

% minimun number of pulses to do further analysis >= p.minPulseNum
fprintf('Pulses to delete: ')
pul2rem = zeros(size(pInf.wc, 1), 1);
if numel(pInf.wc) < p.minPulseNum; fprintf(' 1 (single) '); pul2rem(:) = 1; end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% doing relative amplitude & pLik & distance & stim start thresholding
% stim start in ms vDat.silencePre (prior no pulses can happen prior stimulation)
% FIX: remove - this could be set in paramters up on call??
pul2rem(pInf.pLik < p.pLik_Thrs | pInf.ampR < p.ampR_Thrs | pInf.proxSng > p.p2p_disThrs) = 1;
if sum(pul2rem); fprintf([num2str(sum(pul2rem)),' (low pLik) or (low ampR) or (high dist) ']); end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% minimun number of pulses to do further analysis >= p.minPulseNum
if numel(pInf.wc) < p.minPulseNum; fprintf(' 1 (single) '); pul2rem(:) = 1; end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
[~, ~, pInf] = findSongBouts(pInf, wInf);


% ************************************************************
% ********************** Generate Bouts **********************
% ************************************************************

% delete pulses outside bouts / and small sine bouts using "p.sine_LThrs", "p.minPn", "p.minSl".
[~, boutMask] = findSongBouts(pInf, wInf);
pul2rem(boutMask(pInf.wc) == 0) = 1;
if sum(pul2rem); fprintf([num2str(sum(pul2rem)),' (outside bouts) ']); end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% minimun number of pulses to do further analysis >= p.minPulseNum
if numel(pInf.wc) < p.minPulseNum; fprintf(' 1 (single) '); pul2rem(:) = 1; end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
[~, ~, pInf] = findSongBouts(pInf, wInf);

% ******************************************************************************
% ********************** Grow bouts based on dist and pLik *********************
% ******************************************************************************

% basically add pulses nearby bouts (>= 2 pulses || pLik > 10), uses (p.p2p_growbout_disThrs)
pInf = growBouts(pInf, wInf, dInf);

% minimun number of pulses to do further analysis >= p.minPulseNum
if numel(pInf.wc) < p.minPulseNum; fprintf(' 1 (single) '); pul2rem(:) = 1; end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
[~, ~, pInf] = findSongBouts(pInf, wInf);

% *************************************************************************
% ********************** Eliminate bad single pulses **********************
% *************************************************************************

% pLik and ampR thresholding for single pulses
pul2rem((pInf.pLik < p.pLik_Thrs_sPul | pInf.ampR < p.ampR_Thrs_sPul | pInf.proxSng > p.p2p_disThrs ) & pInf.ppBout == 1) = 1;
if sum(pul2rem); fprintf([num2str(sum(pul2rem)),' (low pLik or ampR for single) ']); end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% minimun number of pulses to do further analysis >= p.minPulseNum
if numel(pInf.wc) < p.minPulseNum; fprintf(' 1 (single) '); pul2rem(:) = 1; end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
[~, ~, pInf] = findSongBouts(pInf, wInf);

% eliminate single pulses dist thresholding for single pulses
pul2rem((pInf.proxSng > p.p2p_disThrs ) & pInf.ppBout == 1) = 1;
if sum(pul2rem); fprintf([num2str(sum(pul2rem)),' (high dist for single) ']); end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% minimun number of pulses to do further analysis >= p.minPulseNum
if numel(pInf.wc) < p.minPulseNum; fprintf(' 1 (single) '); pul2rem(:) = 1; end
[pInf, pul2rem] = remBadPul(pInf, pul2rem);

% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
[~, ~, pInf] = findSongBouts(pInf, wInf);
fprintf(' \n')

% *******************************************************
% ********************** Removing bad sine **************
% *******************************************************

% evaluate cases where sine and pulse coincide
% remove sine in between pulses
% restore sine of pulse is not a true one

% prior on sine size? can not be louder than the loudest pulse?

% fuse sine bouts if they are <=150ms apart and no pulse in between
%tmpSong = g.Song;
%g.nThr = g.nAmp*12;
%sRng = stEn(i,1):stEn(i,2);
%sTest = mean(abs(tmpSong(sRng))>g.nAmp);
%mTest = mean(mean(abs(tmpSong(sRng)))./std(abs(tmpSong(sRng))));
%sTest<g.nThr || mTest < 1;

% fuses sine bouts that are 150 ms apart.
wInf = combineSineSong(pInf, wInf, 1);

% passing bSin to aSin (so basically the number of mics who detected  sine)
wInf = aSin2bSin(pInf, wInf);

% ******************************************************
% ********************** Wrapping **********************
% ******************************************************

% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
[~, ~, pInf] = findSongBouts(pInf, wInf);

% generate bout structure and remove pulses outside bouts (bInf.stEn, bInf.Mask)
[bInf.stEn, bInf.Mask, pInf] = findSongBouts(pInf, wInf, 1);
pul2rem(bInf.Mask(pInf.wc) == 0) = 1;
pInf = remBadPul(pInf, pul2rem);

% final wInf update, use bMsk to clean up sine vector and stEn variable
% update proxSng & wInf & temporal bout fields
pInf = findSongNearPulse(pInf, wInf);
wInf = combineSineSong(pInf, wInf, 1, bInf.Mask); % fuse and filter (using mask)
[bInf.stEn, bInf.Mask, pInf] = findSongBouts(pInf, wInf, 1);

% generates extra bout information (bInf.Bout, bInf.bRng, bInf.Type)
bInf = finalizeBoutInfo(pInf, wInf, bInf);
Song = pInf.Song;%int16(pInf.Song*1000);
% pInf.pSec = int16(pInf.pSec*1000);
pInf = rmfield(pInf, {'Song'; 'exSong'; 'boutTag'; 'ppBout'; 'boutLikMax'; 'boutLikMin'; 'stEn'});
wInf = rmfield(wInf, {'nAmp'; 'bSin'});
% storing pulses prior to filtering for offline evaluation
% REMOVE?
dInf = rmfield(dInf, {'Song'; 'exSong'; 'boutTag'; 'ppBout'; 'boutLikMax'; 'boutLikMin'; 'stEn'});
% dInf.pSec = int16(dInf.pSec*1000);

fprintf(['Final number of : pulses ', num2str(length(pInf.wc)), ' , sine epocs ', num2str(size(wInf.stEn, 1)),...
   ' , bouts found ', num2str(size(bInf.stEn, 1)), ' \n\n'])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [detM, wc, c2SR] = findUniPul(sInf, oneSong)
% basically collects pInf.detM, pInf.wc, pInf.c2SR from sInf, and removes pulse centers
% from edges (+-25 ms) and centers that sample the same pulse based on temporal
% distance between them (dPnts), centers with distances lower than 5 ms are deleted and
% distances within 5:12 ms are evaluated to find the center
sampLength = length(oneSong);
fieldName = 'pulseInfo2'; % can only use pulseInfo2 since `pcndInfo` and `pulseInfo` miss `c2SR`
tempWC = cell(length(sInf),1);
tempc2SR = cell(length(sInf),1);
for i = 1:length(sInf)
   if isfield(sInf(i).(fieldName), 'scmx') % removing scmx
      sInf(i).(fieldName) = rmfield(sInf(i).(fieldName), 'scmx');
   end
   if ~isfield(sInf(i).(fieldName), 'wc')
      tempWC{i,1} = []; tempc2SR{i,1} = []; continue;
   end
   tempWC{i,1} = sInf(i).(fieldName).wc(sInf(i).(fieldName).wc>0);
   if isfield(sInf(i).(fieldName), 'c2SR')
      tempc2SR{i,1} = sInf(i).(fieldName).c2SR(sInf(i).(fieldName).wc>0);
   else
      tempc2SR{i,1} = [];
   end
end

EmptyTrace = zeros(numel(sInf), 1);
for i = 1:length(EmptyTrace)
   EmptyTrace(i, 1) = ~isempty(tempWC{i});
end

if sum(EmptyTrace) > 1 || numel(tempWC) == 1 && sum(EmptyTrace) == 1
   pWid = 251; % time ths (+-25 ms)
   aPul = []; % temporal pulse idx, wc and c2SR
   for i = 1:length(sInf)
      aPul = [aPul;[i*ones(length(tempWC{i}),1), tempWC{i}', tempc2SR{i}']];
   end
   [wcSrt mIdx] = sort(aPul(:,2));
   % discarding pulses at the begining and end (<> time ths 125.5)
   mIdx(wcSrt<ceil(pWid/2)) = [];
   wcSrt(wcSrt<ceil(pWid/2)) = [];
   mIdx(wcSrt>(sampLength-ceil(pWid/2))) = [];
   wcSrt(wcSrt>(sampLength-ceil(pWid/2))) = [];
   mIdx = aPul(mIdx,[1 3]);
   dPnts = diff([0;wcSrt]);
   % step to posibly modify for vPR6
   sPnts = find(dPnts>50); % centers closer than 5 ms are fused (mean time point)
   sPntConflict = find(dPnts>50 & dPnts<125,2); % evaluate center distance > 5 ms && < 12.5 ms
   while ~isempty(sPntConflict)
      % iterative removal of pulse center's that belong to a single pulse (5-12ms apart)
      probPnt = sPntConflict(1);
      sPntIdx = find(sPnts == probPnt);
      % get idx distance to pre & post
      numPrePul = sPnts(sPntIdx) - sPnts(sPntIdx-1);
      if sPntIdx < length(sPnts)
         numPostPul = sPnts(sPntIdx+1) - sPnts(sPntIdx);
      else
         numPostPul = 100000;
      end
      % get max amp from current pulse to pre & post (+-1ms)
      preMax = max(abs(oneSong(wcSrt(sPnts(sPntIdx-1))-10:...
         wcSrt(sPnts(sPntIdx)-1)+10)));
      if sPntIdx < length(sPnts)
         postMax = max(abs(oneSong(wcSrt(sPnts(sPntIdx))-10:...
            wcSrt(sPnts(sPntIdx+1)-1)+10)));
      else
         postMax = 0;
      end
      if length(sPntConflict) == 2 && sPntConflict(1)+1 == sPntConflict(2)
         wcSrt(probPnt) = []; mIdx(probPnt,:) = [];
      elseif numPrePul > numPostPul || ...
            (numPrePul == numPostPul && preMax > postMax)
         % assumes wc(probPnt-1) is the true peak, so it deletes
         % wc(probPnt)and following ones depending on numPostPul
         wcSrt(probPnt:probPnt+numPostPul-1) = [];
         mIdx(probPnt:probPnt+numPostPul-1,:) = [];
      else
         % length(sPntConflict) ~= 2 && numPrePul < numPostPul || preMax < postMax
         % assumes wc(probPnt) is the true peak so it deletes
         % wc(probPnt-1)and preceding ones depending on numPostPul
         wcSrt(probPnt-numPrePul:probPnt-1) = [];
         mIdx(probPnt-numPrePul:probPnt-1,:) = [];
      end
      dPnts = diff([0;wcSrt]);
      sPnts = find(dPnts>50);
      sPntConflict = find(dPnts>50 & dPnts<125,2);
   end
   wc = zeros(length(sPnts),1);
   detM = zeros(length(sPnts),length(sInf));
   c2SR = zeros(length(sPnts),length(sInf));
   % for cases where a 5ms distance is found the wc is calculated as the
   % mean of their individual wc (seems that it happens mostly between mic)
   for i = 1:length(sPnts)
      if i == length(sPnts);
         wc(i,1) = mean(wcSrt(sPnts(i):end));
         detM(i,mIdx(sPnts(i):end,1)) = 1;
         c2SR(i,mIdx(sPnts(i):end,1)) = mIdx(sPnts(i):end,2);
      else
         wc(i,1) = mean(wcSrt(sPnts(i):sPnts(i+1)-1));
         detM(i,mIdx(sPnts(i):sPnts(i+1) - 1,1)) = 1;
         c2SR(i,mIdx(sPnts(i):sPnts(i+1)-1,1)) = mIdx(sPnts(i):sPnts(i+1)-1,2);
      end
   end
   wc = round(wc);
else
   detM = zeros(0,2); wc = zeros(0,1); c2SR = zeros(0,2);
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pInf = rmpulseclose2edge(pInf)
% delete pulses close to trial edges using "p.TimeThrs"
global p%p.TimeThrs
pul2rm = pInf.wc < p.TimeThrs | pInf.wc > (numel(pInf.Song) - p.TimeThrs);
if sum(pul2rm); fprintf([num2str(sum(pul2rm)),' pulses close to edges, ']); end
pInf.detM(pul2rm>0,:) = [];
pInf.wc(pul2rm>0,:) = [];
pInf.c2SR(pul2rm>0,:) = [];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wInf = wInfGen(sInf, pInf, wInf)
% generates wInf, threshold sine by amplitude wInf.aSin(bTest' ~= 0) = 1,
% % delete epocs of sine shorter than p.sine_LThrs (default 70 ms)
global p%p.sine_nThrs p.sine_LThrs
% p.sine_nThrs = 0.15;
oneSong = pInf.Song;
%% generating bSin / aSin
% get sine epocs from both mics
aSin = zeros(size(oneSong)); bSin = zeros(size(oneSong));
wInf.aSin = zeros(size(oneSong)); wInf.bSin = zeros(size(oneSong));
for i = 1:length(sInf)
   if isempty(sInf(i).winSine.NWK3060)
      sInf(i).winSine.NWK3060.num_events = [];
      sInf(i).winSine.NWK3060.start = [];
      sInf(i).winSine.NWK3060.stop = [];
      sInf(i).winSine.NWK3060.length = [];
   end
   if isempty(sInf(i).winSine.NWK2012)
      sInf(i).winSine.NWK2012.num_events = [];
      sInf(i).winSine.NWK2012.start = [];
      sInf(i).winSine.NWK2012.stop = [];
      sInf(i).winSine.NWK2012.length = [];
   end
   
   for j = 1:size(sInf(i).winSine.NWK3060.start,1);
      sRng = round(sInf(i).winSine.NWK3060.start(j):...
         sInf(i).winSine.NWK3060.stop(j));
      aSin(sRng) = aSin(sRng)+1;
   end
   
   for j = 1:size(sInf(i).winSine.NWK2012.start,1);
      sRng = round(sInf(i).winSine.NWK2012.start(j):...
         sInf(i).winSine.NWK2012.stop(j));
      bSin(sRng) = bSin(sRng)+1;
   end
end
% Smooth version of amplitude for sine epoc
bTest = smooth(double(abs(oneSong) > wInf.nAmp).*(bSin > 0),1001) > p.sine_nThrs;
% In this setup we should not necessary look for consistency between mics,
% instead the amplitude ths should have a greater impact on filtering
% out sine-like noise
% wInf.aSin = 2*aSin + bTest'; % not sure why multipying it by 2?
wInf.aSin(bTest' ~= 0) = 1;
% bSin stores # of mics detecting sine
wInf.bSin = bSin;
%% generating stEn and thresholding sine epoc length
tDat = [false, wInf.aSin~=0, false];
stEn = [strfind(tDat, [0,1]); strfind(tDat, [1,0])-1]'; % get sine start / end
dIdx = [];
for i = 1:size(stEn, 1) % evaluate each sine epoc for:
   % delete epocs of sine shorter than p.sine_LThrs (default 70 ms)
   if diff(stEn(i, :)) < p.sine_LThrs;
      dIdx = [i; dIdx]; wInf.aSin(stEn(i, 1):stEn(i, 2)) = 0;
   end
end
stEn(dIdx, :) = [];
if isempty(stEn); stEn = zeros([], 2); end
wInf.stEn = stEn;
wInf.aSin = cast(wInf.aSin, 'int8');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wInf = aSin2bSin(pInf, wInf)
% pass mic number (bSin) to binary aSin
global p%p.sine_nThrs
for i = 1:size(wInf.stEn, 1)
   if sum(abs(pInf.Song(wInf.stEn(i, 1):wInf.stEn(i, 2))) > p.sine_nThrs) > 50
      wInf.aSin(wInf.stEn(i, 1):wInf.stEn(i, 2)) = ...
         wInf.bSin(wInf.stEn(i, 1):wInf.stEn(i, 2));
   end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wInf = combineSineSong(pInf, wInf, blkP, bMsk)
% Always delete epocs of sine shorter than p.sine_LThrs (default 70 ms)
% fuse sine bouts that are <= 150ms apart,
% temporal variables
global p%p.sine_LThrs
tDat = [false, wInf.aSin~=0, false];
stEn = [strfind(tDat, [0,1]); strfind(tDat, [1,0])-1]'; % get sine start / end
%% fuse sine epocs 150 ms apart if no pulses inbetween
i = 1;
while i ~= size(stEn, 1) && size(stEn, 1) ~= 0 && exist('blkP', 'var') &&  blkP == 1;
   if any(pInf.wc < stEn(i+1, 1) & pInf.wc > stEn(i, 2))
      i = i+1; continue;
   end
   if stEn(i+1, 1)-stEn(i, 2) < 1500
      wInf.aSin(stEn(i, 2):stEn(i+1, 1)) = 1;
      stEn(i, 2) = stEn(i+1, 2); stEn(i+1, :) = [];
   else
      i = i+1;
   end
end
%% Delete Sine epocs shorter than 70 ms
dIdx = [];
for i = 1:size(stEn, 1) % evaluate each sine epoc for:
   % delete epocs of sine shorter than p.sine_LThrs (default 70 ms)
   if diff(stEn(i, :)) < p.sine_LThrs;
      dIdx = [i; dIdx]; wInf.aSin(stEn(i, 1):stEn(i, 2)) = 0;
   end
end
stEn(dIdx, :) = [];
%% Use bMsk to clean up sine vector and stEn variable
if exist('bMsk', 'var')
   wInf.aSin(bMsk == 0) = 0;
   for i = 1:length(pInf.wc)
      wInf.aSin(pInf.wc(i) + [-50:50]) = 0;
   end
   tDat = [false, wInf.aSin~=0, false];
   stEn = [strfind(tDat, [0, 1]); strfind(tDat, [1, 0]) - 1]';
end
%% Delete Sine epocs shorter than 70 ms
dIdx = [];
for i = 1:size(stEn, 1) % evaluate each sine epoc for:
   % delete epocs of sine shorter than p.sine_LThrs (default 70 ms)
   if diff(stEn(i, :)) < p.sine_LThrs;
      dIdx = [i; dIdx]; wInf.aSin(stEn(i, 1):stEn(i, 2)) = 0;
   end
end
stEn(dIdx, :) = [];
%% Update output variables
if isempty(stEn); stEn = zeros([], 2); end
wInf.stEn = stEn;
wInf.aSin = cast(wInf.aSin, 'int8');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pInf = fcmx_calculator(pInf, sInf, wInf, mNoiseSamp, noiseOffset, data)
% updates wc based on peak frequency (-+ 7.5 ms from original wc), generates: pInf.dog and pInf.fcmx
% get pSec (pulse trace per mic)
pSec = zeros(251, length(sInf), length(pInf.wc)); % % pSec(-/+12.5ms,mic#,pulse#)
rng = zeros(1, length(pInf.wc)*251); % timepoints +/- 12.5 ms
for i = 1:length(pInf.wc)
   rng(1,(i-1)*251+1:i*251) = pInf.wc(i)-125:pInf.wc(i)+125;
end
intOff = zeros(length(sInf),1);
for i = 1:length(sInf)
   % loading trace per mic
   tempDat = data(i, rng);
   % getting pulse +/- 12.5ms
   tempDat = (double(tempDat) - noiseOffset(i))/1000; % substract noise offset (mean)
   tempDat = reshape(tempDat, 251, 1, length(pInf.wc));
   pSec(:, i, :) = single(tempDat);
   clear tempDat
   % getting whitenning noise and get the mean of abolute deviation
   intOff(i, 1) = mean(abs((double(mNoiseSamp(:, i)) - noiseOffset(i))/1000));
   % integral of a pulse
   if min(size(pSec(:, i, :))) == 0
      pInf.mAmp.Int = zeros(0, length(sInf));
   else
      pInf.mAmp.Int(:, i) = squeeze(trapz((abs(pSec(:, i, :)) - intOff(i)).^2))';
   end
end
for i = 1:length(pInf.wc)
   tempval = squeeze(max(abs(pSec(:, :, i))))';
   if size(tempval, 2) < 2; tempval = tempval'; end
   pInf.mAmp.Mx(i, :) = tempval; % pulse max
   clear tempval
end
if ~isfield(pInf.mAmp, 'Mx'); pInf.mAmp.Mx = zeros(length(pInf.wc), 2); end
pInf.ampR = max(pInf.mAmp.Mx, [], 2)/wInf.nAmp; % wInf.nAmp == nAmp
[~, mIdx] = max(pInf.mAmp.Mx, [],  2);
pSecTemp = zeros(size(pSec, 1), size(pSec, 3));
for i = 1:size(pSec, 3)
   pSecTemp(:,i) = pSec(:,mIdx(i), i);
end
pInf.pSec = squeeze(pSecTemp)';
pInf = rmfield(pInf,{'mAmp','ampR'});
% updating each wc and getting wc
FetchParams % loading params
% assign params
ngw = numel(Params.DoGwvlt); fc = Params.fc ; fs = Params.Fs;
wvlt = cell(1,ngw); width = ceil(size(pInf.pSec,2)/2);
for i = 1:ngw; wvlt{i} = ['gaus' num2str(Params.DoGwvlt(i))]; end
sc = zeros(ngw,numel(fc));
for i = 1:numel(wvlt); sc(i,:) = scales_for_freqs(fc,1/fs,wvlt{i}); end
if ~isempty(pInf.wc)
   for pIdx = 1:numel(pInf.wc)
      fprintf('*')
      % xs == pulse shape
      xs = pInf.pSec(pIdx,:);
      cmhSong = single(zeros(1,numel(xs))); % Storage for the maximum mexican hat wavelet coefficient for each bin.
      cmh_dog = int8(zeros(1,length(xs))); % Storage for the order of the D.o.G. wavelet for which the highest coefficient occured.
      cmh_sc = int8(zeros(1,length(xs))); % Storage for the scale at which the highest mexican hat coefficient occured.
      for i = 1:numel(wvlt)
         for j = 1:size(sc,2)
            temp = single(abs(cwt(xs,sc(i,j),wvlt{i})));
            cmh_sc(temp>cmhSong) = j; cmh_dog(temp>cmhSong) = i;
            cmhSong(temp>cmhSong) = temp(temp>cmhSong);
            clear temp;
         end
      end
      % getting max per pulse
      sig4Test = runningExtreme(cmhSong,Params.pWid,'max');
      sig4Test = smooth(sig4Test,(fs/1000)+Params.pWid);
      % update pInf.wc, using wc (relative idx) only within -+7.5 ms of current center
      wc = find(sig4Test == max(sig4Test(50:200)));
      wc = round(mean(wc)); wc = wc - width;
      pInf.wc(pIdx, 1) = pInf.wc(pIdx) + wc;
      % update dog_at_max fc_at_max sc_at_max
      peak = width + wc;
      dog_at_max = cmh_dog(peak);
      fc_at_max = fc(cmh_sc(peak));
      pInf.dog(pIdx,1) = dog_at_max;
      pInf.fcmx(pIdx,1) = fc_at_max;
      clear xs cmhSong cmh_dog cmh_sc peak dog_at_max
      clear fc_at_max sc_at_max wc peak sig4Test
   end
   fprintf('\n')
else
   pInf.dog = zeros([],1); pInf.fcmx = zeros([],1);
end
pInf = rmfield(pInf,{'pSec'});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pInf = getpsec(pInf, sInf, wInf, noiseOffset, mNoiseSamp, data)
% generates pInf.pSec, amplitude params (pInf.mAmp.Int, pInf.mAmp.Mx, pInf.ampR) and updates pInf.c2SR
% At the end it replace pulses onto song trace (pulse coming from the mic with high amplitude)
% getting amplitude per pulse per mic
pSec = zeros(251, length(sInf), length(pInf.wc)); % % pSec(-/+12.5ms,mic#,pulse#)
rng = zeros(1, length(pInf.wc)*251); % timepoints +/- 12.5 ms
for i = 1:length(pInf.wc)
   rng(1,(i-1)*251+1:i*251) = pInf.wc(i)-125:pInf.wc(i)+125;
end
intOff = zeros(length(sInf), 1);
for i = 1:length(sInf)
   % loading trace per mic
   tempDat = data(i,rng);
   % getting pulse +/- 12.5ms
   tempDat = (double(tempDat) - noiseOffset(i))/1000; % substract noise offset (mean)
   tempDat = reshape(tempDat, 251, 1, length(pInf.wc));
   pSec(:, i, :) = single(tempDat);
   clear tempDat
   % getting whitenning noise and get the mean of abolute deviation
   intOff(i, 1) = mean(abs((double(mNoiseSamp(:, i)) - noiseOffset(i))/1000));
   % integral of a pulse
   if min(size(pSec(:, i, :))) == 0
      pInf.mAmp.Int = zeros(0, length(sInf));
   else
      pInf.mAmp.Int(:, i) = squeeze(trapz((abs(pSec(:, i, :)) - intOff(i)).^2))';
   end
end
for i = 1:length(pInf.wc)
   tempval = squeeze(max(abs(pSec(76:176, :, i))))';
   if size(tempval, 2) < 2; tempval = tempval'; end
   pInf.mAmp.Mx(i, :) = tempval; % pulse max
   clear tempval
end
if ~isfield(pInf.mAmp, 'Mx'); pInf.mAmp.Mx = zeros(length(pInf.wc), 2); end
pInf.ampR = max(pInf.mAmp.Mx, [], 2)/wInf.nAmp; % wInf.nAmp == nAmp
% getting c2SR per pulse
[~, mIdx] = max(pInf.mAmp.Mx, [],  2);
% disregarding instances of abnormal high signal (probably flying)
pSecTemp = zeros(size(pSec, 1), size(pSec, 3));
for i = 1:size(pSec, 3);
   pSecTemp(:,i) = pSec(:, mIdx(i), i);
end
pInf.pSec = squeeze(pSecTemp)';
end

function [pInf] = centerPulse(pInf)
% hard coded to 251 bin pSec
fnorm = @(x) x/norm(x);
if ~isempty(pInf.wc)
   try
      %% normalizing pulses by its vectorial norm
      clear rawDat nDat pLab;
      rawDat = double(pInf.pSec);
      nSamp = size(rawDat, 1);
      for j=1:nSamp
         rawDat(j, :) = fnorm(rawDat(j, :)); %normalize to unity
      end
      %% Center on max power and flip if negative lobe to left of max power
      oLen = length(rawDat(1, :));
      nLen = 2*oLen;
      nDat = zeros(nSamp, nLen);
      for j = 1:nSamp
         %% centering to max power
         tPul = rawDat(j, :);
         [~, mIdx(j, 1)] = max(smooth(tPul.^2, 15)); %get max power, then center on this
         lPad = oLen - mIdx(j, 1);  %i.e. ((total_length/2) - C)
         rPad = nLen - oLen - lPad;
         nDat(j, :) = [zeros(1, lPad)  pInf.pSec(j, :)  zeros(1, rPad)];
         %% flipping sign
         if mean(nDat(j, oLen-10:oLen)) <0; nDat(j, :) = nDat(j, :)*-1; end
         %% correct wc
         pInf.wc(j, 1) = pInf.wc(j, 1) - (mIdx(j, 1) - 125) + 1;
      end
      nDat = nDat(:, 126:376);
      pInf.pSec = nDat;
      pul2rem = mIdx < 76 | mIdx > 176; % delete all pulses that do not have a clear center (only allows +-5ms)
      [pInf, ~] = remBadPul(pInf, pul2rem);
   catch ME
      disp(ME.getReport())
      nDat = [];% why??
   end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf, wInf] = getp_Amp_c2SR(pInf, sInf, wInf, noiseOffset, mNoiseSamp, data)
% generates pInf.pSec, amplitude params (pInf.mAmp.Int, pInf.mAmp.Mx, pInf.ampR) and updates pInf.c2SR
% At the end it replaces pulses in song trace with the pulse coming from the mic with high amplitude
% getting amplitude per pulse per mic
pSec = zeros(251, length(sInf), length(pInf.wc)); % % pSec(-/+12.5ms,mic#,pulse#)
rng = zeros(1, length(pInf.wc)*251); % timepoints +/- 12.5 ms
for i = 1:length(pInf.wc)
   rng(1,(i-1)*251+1:i*251) = pInf.wc(i)-125:pInf.wc(i)+125;
end
intOff = zeros(length(sInf), 1);
for i = 1:length(sInf)
   tempDat = data(i,rng);
   % getting pulse +/- 12.5ms
   tempDat = (double(tempDat) - noiseOffset(i))/1000; % substract noise offset (mean)
   tempDat = reshape(tempDat, 251, 1, length(pInf.wc));
   pSec(:, i, :) = single(tempDat);
   % storing mic data (m*)
   clear tempDat
   % getting whitenning noise and get the mean of abolute deviation
   intOff(i, 1) = mean(abs((double(mNoiseSamp(:, i)) - noiseOffset(i))/1000));
   % integral of a pulse
   if min(size(pSec(:, i, :))) == 0
      pInf.mAmp.Int = zeros(0, length(sInf));
   else
      pInf.mAmp.Int(:, i) = squeeze(trapz((abs(pSec(:, i, :)) - intOff(i)).^2))';
   end
end
for i = 1:length(pInf.wc)
   tempval = squeeze(max(abs(pSec(76:176, :, i))))';
   if size(tempval, 2) < 2; tempval = tempval'; end
   pInf.mAmp.Mx(i, :) = tempval; % pulse max
   clear tempval
end
if ~isfield(pInf.mAmp, 'Mx'); pInf.mAmp.Mx = zeros(length(pInf.wc), 2); end
pInf.ampR = max(pInf.mAmp.Mx, [], 2)/wInf.nAmp; % wInf.nAmp == nAmp
% getting c2SR per pulse
[~, mIdx] = max(pInf.mAmp.Mx, [],  2);
temp = [];
for i = 1:length(mIdx); temp(i,1) = pInf.c2SR(i,mIdx(i)); end
pInf.c2SR = temp; clear temp
if isempty(pInf.c2SR); pInf.c2SR = zeros(length(pInf.wc), 1); end
% disregarding instances of abnormal high signal (probably flying)
pSecTemp = zeros(size(pSec, 1), size(pSec, 3));
for i = 1:size(pSec, 3);
   pSecTemp(:,i) = pSec(:,mIdx(i), i);
   % update c2SR for all
   sidAmps =[max(abs(pSecTemp(1:50, i))) max(abs(pSecTemp(201:251, i)))];
   pInf.c2SR(i, 1) = max(abs(pSecTemp(90:160, i)))/max(sidAmps);
end
pInf.pSec = squeeze(pSecTemp)';
% replacing pulse onto song trace (pulse coming from the mic with high amplitude)
for i = 1:length(pInf.wc);
   pInf.Song(pInf.wc(i) + [-125:125]) = pInf.pSec(i,:);
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pInf = getp_hpWi_pulse2song(pInf)
pSec = pInf.pSec(:, 51:201); % pSec(-/+7.5ms,pulse#) 150 timepoints
pInf.hpWi = zeros(size(pSec, 1), 1);
for i = 1:size(pSec, 1)
   testPul = resample(pSec(i, :), 10, 1);
   corFact = sign(testPul(abs(testPul) == max(abs(testPul))));
   testPul = smooth(testPul,100);
   if length(corFact) > 1 % this happens rarely
      testPul = testPul*sign(corFact(1));
      pInf.hpWi(i, 1) = getHPWPerPulse(testPul);
   else
      testPul = testPul*sign(corFact);
      try
         pInf.hpWi(i, 1) = getHPWPerPulse(testPul);
      catch ME
         disp(ME.getReport)
         pInf.hpWi(i,1) = length(testPul);
      end
   end
end
end

function [hpWi] = getHPWPerPulse(pulse)
if ~isrow(pulse); pulse = pulse'; end
pulse(pulse < max(pulse)./sqrt(2)) = 0;
pulse(pulse > max(pulse)./sqrt(2)) = 1;
hpWi = diff([0 pulse]);
hpWi = diff([find(hpWi == 1,1) find(hpWi == -1,1, 'last')]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf, pul2rem] = remBadPul(pInf, pul2rem)
% removes labeled pulses on pul2rem, also re-sort all values based on wc
pInf.wc(pul2rem>0,:) = [];
if isfield(pInf, 'detM'); pInf.detM(pul2rem>0,:) = []; end
if isfield(pInf, 'c2SR'); pInf.c2SR(pul2rem>0,:) = []; end
if isfield(pInf, 'mAmp'); pInf.mAmp.Mx(pul2rem>0,:) = []; pInf.mAmp.Int(pul2rem>0,:) = []; end
if isfield(pInf, 'ampR'); pInf.ampR(pul2rem>0,:) = []; end
if isfield(pInf, 'pSec'); pInf.pSec(pul2rem>0,:) = []; end
if isfield(pInf, 'hpWi'); pInf.hpWi(pul2rem>0,:) = []; end
if isfield(pInf, 'Model'); pInf.Model(pul2rem>0,:) = []; end
if isfield(pInf, 'dog'); pInf.dog(pul2rem>0,:) = []; end
if isfield(pInf, 'fcmx'); pInf.fcmx(pul2rem>0,:) = []; end
if isfield(pInf, 'pLik'); pInf.pLik(pul2rem>0,:) = []; end
if isfield(pInf, 'proxSng'); pInf.proxSng(pul2rem>0,:) = []; end
if isfield(pInf, 'boutTag'); pInf.boutTag(pul2rem>0,:) = []; end
if isfield(pInf, 'ppBout'); pInf.ppBout(pul2rem>0,:) = []; end
pul2rem = 0*pInf.wc;
% sort remaining pulses
if numel(pInf.wc) > 0
   [~, srtIdx] = sort(pInf.wc);
   pInf.wc = pInf.wc(srtIdx,:);
   if isfield(pInf, 'detM'); pInf.detM = pInf.detM(srtIdx,:); end
   if isfield(pInf, 'c2SR'); pInf.c2SR = pInf.c2SR(srtIdx,:); end
   if isfield(pInf, 'mAmp'); pInf.mAmp.Mx = pInf.mAmp.Mx(srtIdx,:); pInf.mAmp.Int = pInf.mAmp.Int(srtIdx,:); end
   if isfield(pInf, 'ampR'); pInf.ampR = pInf.ampR(srtIdx,:); end
   if isfield(pInf, 'pSec'); pInf.pSec = pInf.pSec(srtIdx,:); end
   if isfield(pInf, 'hpWi'); pInf.hpWi = pInf.hpWi(srtIdx,:); end
   if isfield(pInf, 'Model'); pInf.Model = pInf.Model(srtIdx,:); end
   if isfield(pInf, 'dog'); pInf.dog = pInf.dog(srtIdx,:); end
   if isfield(pInf, 'fcmx'); pInf.fcmx = pInf.fcmx(srtIdx,:); end
   if isfield(pInf, 'pLik'); pInf.pLik = pInf.pLik(srtIdx,:); end
   if isfield(pInf, 'proxSng'); pInf.proxSng = pInf.proxSng(srtIdx,:); end
   if isfield(pInf, 'boutTag'); pInf.boutTag = pInf.boutTag(srtIdx,:); end
   if isfield(pInf, 'ppBout'); pInf.ppBout = pInf.ppBout(srtIdx,:); end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf] = findSongNearPulse(pInf, wInf)
% this function finds the minimun distance to other pulses / sine
% requires at least 2 pulses
pInf.proxSng = 0*pInf.wc;
tempSin = find(wInf.aSin)';
if numel(pInf.wc) > 1
   for c_i = 1:numel(pInf.wc)
      if c_i == 1
         tstPul = [pInf.wc(c_i+1); tempSin(1:200:end)];
      elseif c_i == numel(pInf.wc)
         tstPul = [pInf.wc(c_i-1); tempSin(1:200:end)];
      else
         tstPul = [pInf.wc(c_i-1); pInf.wc(c_i+1); tempSin(1:200:end)];
      end
      tPul = abs(tstPul - pInf.wc(c_i));
      pInf.proxSng(c_i,1) = min(tPul);
   end
elseif numel(pInf.wc) == 1
   pInf.proxSng(1,1) = Inf;
else
   pInf.proxSng = zeros([],1);
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stEn, boutMask, pInf] = findSongBouts(pInf, wInf, fin)
global p%p.TimeThrs p.s2s_disThrs p.minPn p.minSl
%% Initialize sBouts
sBouts = pInf.Song*0; sBouts(pInf.wc) = 1; sBouts(wInf.aSin > 0) = 1;
% disregard signal from 0:p.TimeThrs or end-p.TimeThrs:end
sBouts([1:(p.TimeThrs+1), (end-p.TimeThrs):end]) = 0;
tDat = [false, sBouts ~= 0, false];
stEn = [strfind(tDat, [0,1]); strfind(tDat, [1,0])-1]';
sBouts(pInf.wc) = 10;
sBouts(wInf.aSin > 0) = 2;
%% Bout definition: as pulses/sine closer than "p.s2s_disThrs" (us)
idx = 1;
while idx < size(stEn,1);
   if stEn(idx+1,1)-stEn(idx,2) < p.s2s_disThrs
      sBouts(stEn(idx,2)+1:stEn(idx+1,1)-1) = 1;
      stEn(idx,:) = [stEn(idx,1) stEn(idx+1,2)];
      stEn(idx+1,:) = [];
   else
      idx = idx + 1;
   end
end
%% delete bouts shorter than < 70 ms, if it has sine on it.
if ~isempty(stEn)
   for i = 1:size(stEn,1)
      sine_in(i,1) = sum(sBouts(stEn(i,1):stEn(i,2)) == 2) > 0;
   end
   stEn((diff(stEn, [], 2) < 700 & sine_in),:) = [];
   clear sine_in
end
%% delete bouts close to end/start using p.TimeThrs / low pulse number (3,2) and short sine (less than 300 ms)
dIdx = [];
for i = 1:size(stEn,1)
   bRng = stEn(i,1):stEn(i,2);               % bout mask for current bout
   nPulses = sum(sBouts(bRng) == 10);        % number of pulses for current  bout
   nSine = (sum(sBouts(bRng) == 2))/1000;    % ms of sine song in current bout
   snRng = pInf.Song(bRng);                  % song trace for current bout
   if bRng(1) < p.TimeThrs || bRng(end) > length(pInf.Song)-p.TimeThrs || ...
         (nPulses < p.minPn(1) && nSine < p.minSl(1)) || (nPulses < p.minPn(2) && nSine < p.minSl(2))
      dIdx = [dIdx ; i];
   elseif exist('fin', 'var') && ((nPulses < p.minPn(3) && nSine < p.minSl(3)) || ...
         (nPulses == p.minPn(4) && nSine < p.minSl(4)))
      dIdx = [dIdx ; i];
   elseif nPulses < p.minPn(6) && exist('fin', 'var')
      sTest = mean(abs(snRng(sBouts(bRng) == 2)) > wInf.nAmp);
      % at least 10% of sine should be > wInf.nAmp
      if (nPulses < p.minPn(5) && nSine < p.minSl(5)) || sTest < 0.1
         dIdx = [dIdx ; i];
      end
   end
end
stEn(dIdx,:) = [];
%% get boutMask
boutMask = int8(pInf.Song*0 - 1);
for i = 1:size(stEn,1)
   boutMask(stEn(i,1):stEn(i,2)) = 1;
end
boutMask(wInf.aSin > 0) = boutMask(wInf.aSin > 0) + 1;
boutMask(boutMask == -1) = 0;
if exist('fin', 'var')
   idx = 1;
   while idx < size(stEn,1)
      % fusing bouts that are closer than p.s2s_disThrs (default 200 ms)
      if stEn(idx+1,1)-stEn(idx,2) < p.s2s_disThrs
         boutMask(stEn(idx,2)+1:stEn(idx+1,1)-1) = 1;
         stEn(idx,:) = [stEn(idx,1) stEn(idx+1,2)];
         stEn(idx+1,:) = [];
      else
         idx = idx + 1;
      end
   end
end
%% copy bout info on pInf
pInf.boutTag = zeros(numel(pInf.wc), 1);
pInf.ppBout = zeros(numel(pInf.wc), 1);
pInf.boutLikMax = zeros(numel(pInf.wc), 1);
pInf.boutLikMin = zeros(numel(pInf.wc), 1);
if isempty(stEn)
   stEn = zeros([], 2);
end
for bout_i = 1:size(stEn, 1)
   pIdxperBout = pInf.wc <= stEn(bout_i, 2) & pInf.wc >= stEn(bout_i, 1);
   pInf.boutTag(pIdxperBout, 1) = bout_i;
   pInf.ppBout(pIdxperBout, 1) = sum(pIdxperBout);
   pInf.boutLikMax(pIdxperBout, 1) = max(pInf.pLik(pIdxperBout, 1));
   pInf.boutLikMin(pIdxperBout, 1) = min(pInf.pLik(pIdxperBout, 1));
end
pInf.stEn = stEn;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bInf] = finalizeBoutInfo(pInf, wInf, bInf)
global p%p.buff
if isempty(bInf.stEn)
   bInf.stEn = [];
   bInf.Mask = [];
   bInf.p.buff = [];
   bInf.Bout = [];
   bInf.bRng = [];
   bInf.bMsk = [];
   bInf.FrmS = [];
   bInf.bFrm = [];
   bInf.xyrs = [];
   bInf.xyst = [];
   bInf.Type = [];
   bInf.fJmp = [];
   return
end

% SongSeg_getfName;
%p.buff = [35000 35000];
tempSin = find(wInf.aSin > 0)';
temp = 0*bInf.Mask;
temp(pInf.wc) = 1;
temp(tempSin) = 2;
iPnts =  find(temp>0);
if numel(iPnts) ~= 1
   tDat = interp1(iPnts, double(temp(iPnts)),...
      1:length(temp), 'nearest');
else
   tDat(temp == 0) = nan;
end
bInf.Mask = tDat.*double(bInf.Mask>0);
bInf.Mask = int8(bInf.Mask);
Song = pInf.Song;
bInf.p.buff = p.buff;
for i = 1:size(bInf.stEn,1)
   bInf.Bout{i,1} = Song(bInf.stEn(i,1)-p.buff(1):bInf.stEn(i,2)+p.buff(2));
   if max(bInf.Bout{i,1})<3
      bInf.Bout{i,1} = int16(bInf.Bout{i,1}*1000);
   end
   bInf.bRng(i,:) = [p.buff(1)+1 (length(bInf.Bout{i}) - p.buff(2))];
   bInf.bMsk{i,1} = bInf.Mask(bInf.stEn(i,1):bInf.stEn(i,2))>1;
   if isempty(find(bInf.bMsk{i,1} == 1, 1))
      bInf.Type{i,1} = 'Pul';
   elseif isempty(find(bInf.bMsk{i,1} == 0,1))
      bInf.Type{i,1} = 'Sin';
   else
      bInf.Type{i,1} = 'Mix';
   end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf] = growBouts(pInf, wInf, dInf)
% TODO: only does that for pulses - implement for sine as well
fprintf(' growing bouts ')
global p%p.p2p_growbout_disThrs p.s2s_disThrs p.pLik_Thrs_nPul
gate = 1;
pulsesAdded = 0;
while gate == 1
   %% getting pulses that have been rejected
   dInf = dInf_Minus_pInf(dInf, pInf);
   %% determine if there are rejected pulses within the pLik Ths around current bouts
   [bIdx, localIdx] = unique(pInf.boutTag(:, 1));
   pPerBout = pInf.ppBout(localIdx, 1);
   %pLikPerBout = pInf.boutLikMax(localIdx, 1);
   bIdx = bIdx(pPerBout >= 1); % only allow bouts of at least 1 pulses to grow
   gate = 0; % stop if it does not find any nearby pulses
   if ~isempty(bIdx)
      for boutIdx = bIdx'
         nPulse = ((dInf.wc(:, 1) >= (pInf.stEn(boutIdx, 1) - p.p2p_growbout_disThrs) & ...
            dInf.wc(:, 1) <= (pInf.stEn(boutIdx, 2) + p.p2p_growbout_disThrs)) ...
            & dInf.pLik > p.pLik_Thrs_nPul(1)) | ((dInf.wc(:, 1) >= (pInf.stEn(boutIdx, 1) - p.s2s_disThrs) & ...
            dInf.wc(:, 1) <= (pInf.stEn(boutIdx, 2) + p.s2s_disThrs)) ...
            & dInf.pLik >= p.pLik_Thrs_nPul(2));
         nPul = nPulse > 0;
         % add extra pulses to pInf
         if sum(nPul) > 0
            pulsesAdded = pulsesAdded + sum(nPul);
            pInf = pInfAdder(pInf, dInf, nPul);
            dInf = dInf_Minus_pInf(dInf, pInf);
         end
         clear nPulse2 nPulse1 nPul
      end
   end
   pInf = findSongNearPulse(pInf, wInf); % update proxSng
   [~, ~, pInf] = findSongBouts(pInf, wInf); % update temporal bout fields
end
fprintf(' %d pulses added , ', pulsesAdded)
end

function dInf = dInf_Minus_pInf(dInf, pInf)
[~, selPulses] = setdiff(dInf.wc(:, 1), pInf.wc(:, 1));
pul2rem = ones(size(dInf.wc, 1), 1);
pul2rem(selPulses, 1) = 0;
[dInf] = remBadPul(dInf, pul2rem);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pInf = pInfAdder(pInf, dInf, nPul)
% adds pulses stored in dInf to pInf
pInf.wc((end + 1):(end + sum(nPul)), :) = dInf.wc(nPul, :);
pInf.detM((end + 1):(end + sum(nPul)), :) = dInf.detM(nPul, :);
pInf.c2SR((end + 1):(end + sum(nPul)), :) = dInf.c2SR(nPul, :);
pInf.mAmp.Mx((end + 1):(end + sum(nPul)), :) = dInf.mAmp.Mx(nPul, :);
pInf.mAmp.Int((end + 1):(end + sum(nPul)), :) = dInf.mAmp.Int(nPul, :);
pInf.ampR((end + 1):(end + sum(nPul)), :) = dInf.ampR(nPul, :);
pInf.pSec((end + 1):(end + sum(nPul)), :) = dInf.pSec(nPul, :);
pInf.hpWi((end + 1):(end + sum(nPul)), :) = dInf.hpWi(nPul, :);
pInf.Model((end + 1):(end + sum(nPul)), :) = dInf.Model(nPul, :);
if isfield(pInf, 'dog'); pInf.dog((end + 1):(end + sum(nPul)), :) = dInf.dog(nPul, :); end
if isfield(pInf, 'fcmx'); pInf.fcmx((end + 1):(end + sum(nPul)), :) = dInf.fcmx(nPul, :); end
if isfield(pInf, 'pLik'); pInf.pLik((end + 1):(end + sum(nPul)), :) = dInf.pLik(nPul, :); end
if isfield(pInf, 'proxSng'); pInf.proxSng((end + 1):(end + sum(nPul)), :) = dInf.proxSng(nPul, :); end
% sort pulses
if numel(pInf.wc) > 0
   [~, srtIdx] = sort(pInf.wc);
   pInf.wc = pInf.wc(srtIdx,:);
   pInf.detM = pInf.detM(srtIdx,:);
   pInf.c2SR = pInf.c2SR(srtIdx,:);
   pInf.mAmp.Mx = pInf.mAmp.Mx(srtIdx,:);
   pInf.mAmp.Int = pInf.mAmp.Int(srtIdx,:);
   pInf.ampR = pInf.ampR(srtIdx,:);
   pInf.pSec = pInf.pSec(srtIdx,:);
   pInf.hpWi = pInf.hpWi(srtIdx,:);
   if isfield(pInf, 'Model'); pInf.Model = pInf.Model(srtIdx,:); end
   if isfield(pInf, 'dog'); pInf.dog = pInf.dog(srtIdx,:); end
   if isfield(pInf, 'fcmx'); pInf.fcmx = pInf.fcmx(srtIdx,:); end
   if isfield(pInf, 'pLik'); pInf.pLik = pInf.pLik(srtIdx,:); end
   if isfield(pInf, 'proxSng'); pInf.proxSng = pInf.proxSng(srtIdx,:); end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ModelIdx, pLik] = sInfFitMultiPulseMod(pSec, pMod, pLik_scaletype)
% function perform pulse model fitting and provides the model Idx
if ~isempty(pSec)
   pNum = size(pSec, 1);
   for Mod_i = 1:numel(pMod)
      [~, Lik2col] = getPulseLikelihood(pSec, pMod(Mod_i), pLik_scaletype);
      Lik2colAll(:, Mod_i)  = Lik2col.LLR_best(:)';
      ModelIdx(:, Mod_i) = Lik2col.ModIdx + (Mod_i-1)*2;
      clear Lik2col
   end
   %% using the model that fits better
   [pLik, Idx] = max(Lik2colAll, [], 2);
   ModelIdx = ModelIdx(Idx);
else
   pLik = []; ModelIdx = [];
end
end
