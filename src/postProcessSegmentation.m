function	[sInf, pInf, wInf, bInf, Song] = postProcessSegmentation(sInf, data, oneSong, mNoiseSamp, dataScalingFactor, lastSampleToProcess, pulseModelName)
%[sInf, pInf, wInf, bInf, Song] = postProcessSegmentation(sInf, oneSong, mNoiseSamp, data, dataScalingFactor)

% TODO:
% - move all dataScaling and casting to double out of this function

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
end
if ~exist('pulseModelName','var') || isempty(pulseModelName)
   pulseModelName = ''; % default to empty
end

% load pulse models
load pulseModels;
if strfind(pulseModelName, 'NM91')
   pModIdx = cellfun(@(x) strcmp('NM91', x), {pMod.fStr}');
else
   pModIdx = cellfun(@(x) strcmp(pulseModelName, x), {pMod.fStr}');
end
if any(pModIdx)
   pMod = pMod(pModIdx);
else
   warning('did not find pulse model matching "%s". defaulting to NM91.', pulseModelName)
   pMod = pMod(cellfun(@(x) strcmp('NM91', x), {pMod.fStr}'));
end

% estimate noise level
noiseOffset = mean(mNoiseSamp);
nAmp = abs(bsxfun(@minus, double(mNoiseSamp), noiseOffset))/1000;    % remove noise mean
for i = 1:size(nAmp,2)
   nAmp(:,i) = runningExtreme(nAmp(:,i),101,'max');                  % for each channel(?) get sliding window max of noise
end
nAmp = max(mean(nAmp));                                              % max of mean over time(?)
wInf.nAmp = nAmp;

% gather pulse information across all channels - remove duplicates
pInf.Song = oneSong';
[pInf.detM, pInf.wc, pInf.c2SR] = findUniPul(sInf, pInf.Song);

% combine sine song across channels
[wInf] = combineSineSong(sInf, pInf, wInf);
delIdx = find(pInf.wc>lastSampleToProcess | ...
   pInf.wc<125 | pInf.wc>(length(pInf.Song)-125));
pInf.detM(delIdx,:) = [];
pInf.wc(delIdx) = [];
wInf.aSin(ceil(lastSampleToProcess):end) = 0;

% process pulses
pSec = zeros(251, length(sInf), length(pInf.wc));                 % init pulse shapes
rng = zeros(1, length(pInf.wc)*251);                              % list of indices for each pulse shape
for i = 1:length(pInf.wc)
   rng(1,(i-1)*251+1:i*251) = pInf.wc(i)-125:pInf.wc(i)+125;
end
intOff = zeros(length(sInf),1);
% allMDat = [];
for i = 1:length(sInf)
   tempDat = double(data(rng,i))./dataScalingFactor;
   if max(tempDat(:))>3 % ?? I think this discriminates different version by peak amp of recording
      tempDat = (double(tempDat) - noiseOffset(i))/1000;
   else
      tempDat = double(tempDat - noiseOffset(i)/1000);
   end
   intOff(i,1) = mean(abs((double(mNoiseSamp(:,i)) - noiseOffset(i))/1000));
   tempDat = reshape(tempDat, 251, 1, length(pInf.wc));
   pSec(:,i,:) = single(tempDat);
   pInf.mAmp.Int(:,i) = squeeze(trapz((abs(pSec(:,i,:)) - intOff(i)).^2))'; % rms amp of each peak
end
pInf.mAmp.Mx = reshape(max(abs(pSec),[],1), length(sInf), size(pSec,3))'; % peak of abs of each peak (amplitude estimate #2)

[~, mIdx] = max(pInf.mAmp.Mx, [],  2);
for i = 1:length(mIdx)
   temp(i,1) = pInf.c2SR(i,mIdx(i));
end
pInf.c2SR = temp;
pul2rem = zeros(size(pSec, 3),1);

data = double(data(:,1:length(sInf)))./dataScalingFactor;
nDat = sum(abs(data)>max(abs(data(:))*0.3),2);
clear data;
%%
pSecTemp = zeros(size(pSec, 1), size(pSec, 3));
blkP = [];
for i = 1:size(pSec,3)                                                     % mark pulses to remove:
   if sum(nDat([-125:-80 80:125]+pInf.wc(i)))>5 && ...                     % bad SNR (?)
         sum(nDat(pInf.wc(i)-10:pInf.wc(i)+10))
      pul2rem(pInf.wc>=pInf.wc(i)-100 & pInf.wc<pInf.wc(i)+1200) = 1;      % previous pulse within 10ms
      % and following pulse within 120ms(?)
      % but this works since we go over all subsequent pulses
      wInf.aSin(pInf.wc(i)-100:pInf.wc(i)+1200) = 0;
      blkP = [blkP; pInf.wc(i)];
   end
   pSecTemp(:,i) = pSec(:,mIdx(i),i);
   if pInf.c2SR(i) == 0
      sidAmps = [max(abs(pSecTemp(1:50,i))) max(abs(pSecTemp(201:251,i)))];
      pInf.c2SR(i,1) = max(abs(pSecTemp(90:160,i)))/max(sidAmps);
   end
end

blkP = [blkP; length(pInf.Song)];
pSec = squeeze(pSecTemp)';
%%
relTest = pMod.shM*sign(pMod.shM(pMod.shM == max(abs(pMod.shM))));
relVal = getHPWForPulse(resample(relTest, 10, 1));
pInf.pSec = pSec;
pSec = pSec(:,51:201);
for i = 1:size(pSec,1)
   testPul = resample(pSec(i,:),10,1);
   corFact = sign(testPul(abs(testPul)== max(abs(testPul))));
   testPul = testPul*sign(mean(corFact));
   if max(pInf.mAmp.Mx(i,:)) < 4*nAmp || length(corFact) > 1
      pInf.hpWi(i,1) = nan;
   else
      testPul = testPul*sign(testPul(abs(testPul)== max(abs(testPul))));
      pInf.hpWi(i,1) = getHPWForPulse(testPul)/relVal;
   end
end
clear pSec* oneSong time temp* ah bh i j nPnts thresh;
for i = 1:length(pInf.wc)
   pInf.Song(pInf.wc(i)+(-125:125)) = pInf.pSec(i,:);
end
pInf.ampR = max(pInf.mAmp.Mx,[],2)/nAmp;

[pInf,pul2rem] = remBadPul(pInf, pul2rem);
pInf.exSong = runningExtreme(abs(pInf.Song), 65, 'max');


%% ELIMINATE REALLY BAD PULSES
[wInf] = combineSineSong(sInf, pInf,wInf);
[pInf,pul2rem] = remWorstPul(pInf, wInf, pul2rem);
origWInf = wInf;
[pInf, wInf, pul2rem] = remSinConPul(pInf, wInf, sInf, pul2rem);
origPInf = pInf;
[pInf,pul2rem] = remWorstPul(pInf, wInf, pul2rem);
if ~isempty(setdiff(origPInf.wc, pInf.wc))
   nIdx = setdiff(origPInf.wc, pInf.wc);
   for i = length(nIdx)
      idx = origPInf.wc == nIdx(i);
      if origWInf.aSin(nIdx(i))>0 && origPInf.ampR(idx) < 5
         wInf.aSin(nIdx(i) + [-125:125]) = ...
            origWInf.aSin(nIdx(i) + [-125:125]);
      end
   end
end

%% fit pulses to model
if numel(pInf.pSec)>0
   [~,pInf.pLik] = fitPulse2GivenModel(pInf.pSec, pMod);
   pInf.pLik = pInf.pLik.LLR_best;
end
pInfOrig = pInf;
pul2rem(pInf.pLik < -8) = 1; % remove unlikely pulses
[pInf,pul2rem] = remBadPul(pInf, pul2rem);

[~, boutMask] = findSongBouts(pInf, wInf, nAmp, blkP);
pul2rem(boutMask(pInf.wc) == 0) = 1;
[pInf,pul2rem] = remBadPul(pInf, pul2rem);

if numel(pInf.pSec)>0
   [~,pInf.pLik] = fitPulse2GivenModel(pInf.pSec, pMod);
   pInf.pLik = pInf.pLik.LLR_best;
end
pInf = findSongNearPulse(pInf,wInf);
moveOn = 0;
while moveOn == 0
   moveOn = 1;
   if numel(pInf.pLik)<2
      continue
   end
   for i = 2:length(pInf.wc)-1
      [nDis nIdx] = sort(abs(pInf.wc - pInf.wc(i)));
      pSng = pInf.proxSng;
      type = length(setdiff([i-1 i+1], nIdx(2:3)));
      type = 1+(type - 1)*-1;
      nPul = abs([pInf.wc(i+1) pInf.wc(i-1)] - pInf.wc(i));
      nPul = sum(nPul>250 & nPul< 800);
      if (pInf.pLik(i)<3 && (pSng(i)<250 || pSng(i)>800)) ||...
            (pInf.pLik(i)<-1 && (pInf.ampR(i)<10 && nPul<type)) || ...
            (pInf.ampR(i)<5 && nPul == 0) || (pInf.pLik(i) < -4 && nPul<type) || ...
            (pInf.ampR(i)>15 && pInf.pLik(i)< 0 && nPul<type) || ...
            (pInf.pLik(i)<0 && (pInf.c2SR(i)<2.5 || pInf.hpWi(i)>5 ||...
            pInf.hpWi(i)<0.5)) || (pInf.hpWi(i) > 4 && pInf.ampR(i) > 40);
         pul2rem(i) = 1;
         moveOn = 0;
      end
   end
   [pInf,pul2rem] = remBadPul(pInf, pul2rem);
   pInf = findSongNearPulse(pInf,wInf);
end

[pInf,pul2rem] = remBadPul(pInf, pul2rem);
pInf = findSongNearPulse(pInf,wInf);

pInfTmp = pInf;
pul2rem(pInf.pLik < -3 | (pInfTmp.proxSng>800 & pInfTmp.ampR <5) | ...
   pInfTmp.hpWi<0.5 | pInfTmp.hpWi> 4 |  pInfTmp.c2SR < 2)=1;
[pInfTmp] = remBadPul(pInf, pul2rem);
[~,pInfTmp.pLik] = fitPulse2GivenModel(pInfTmp.pSec, pMod);
pInfTmp = findSongNearPulse(pInfTmp,wInf);
pInfTmp.pLik = pInfTmp.pLik.LLR_best;
moveOn = 0;
while moveOn ==0
   moveOn = 1;
   for i = 1:length(pInfTmp.wc)
      [nDis, nIdx] = sort(abs(pInfTmp.wc - pInfTmp.wc(i)));
      if pInfTmp.proxSng(i)>800 || pInfTmp.hpWi(i)>5 ...
            pInfTmp.ampR(nIdx(2))>3*pInfTmp.ampR(i) || ...
            pInfTmp.c2SR(i)<2 || pInfTmp.hpWi(i)<0.5 || ...
            (pInf.c2SR(i)<2.5 && pInf.ampR(i)>15 && pInf.c2SR(i)>1);
         pul2rem(pInf.wc == pInfTmp.wc(i)) = 1;
         moveOn = 0;
      end
   end
   [pInfTmp] = remBadPul(pInf, pul2rem);
   pInfTmp = findSongNearPulse(pInfTmp,wInf);
end
[stEn, boutMask] = findSongBouts(pInfTmp, wInf, nAmp, blkP);

regTest = [];
for i = 1:size(stEn,1)
   regTest = [regTest;[stEn(i,1):750:stEn(i,2)-850]'; stEn(i,2)-650];
end
moveOn = 0;s
while moveOn == 0
   moveOn = 1;
   remPul = find(pul2rem);
   tDat = 0*remPul;
   for i = 1:length(remPul)
      if min(abs(pInf.wc(remPul(i)) - regTest)) < 1500;
         regTest = [regTest; pInf.wc(remPul(i))];
         pul2rem(remPul(i)) = 0;
         moveOn = 0;
      end
   end
end
[pInf,pul2rem] = remBadPul(pInf, pul2rem);
pInf = findSongNearPulse(pInf,wInf);
[~, boutMask] = findSongBouts(pInf, wInf, nAmp, blkP);
pul2rem(boutMask(pInf.wc) == 0) = 1;
[pInf,pul2rem] = remBadPul(pInf, pul2rem);
[wInf] = combineSineSong(sInf, pInf, wInf, blkP, 1);
[stEn, boutMask] = findSongBouts(pInf, wInf, nAmp, blkP, 1);
moveOn = 0;
nPnts = find(nDat);
while moveOn == 0
   moveOn = 1;
   for i = 1:size(stEn,1)
      idx = find(nPnts<stEn(i,2), 1, 'last');
      if boutMask(stEn(i,2))==2 && ~isempty(idx) && ...
            stEn(i,2)-nPnts(idx)<100
         wInf.aSin(nPnts(idx)-100:stEn(i,2)) = 0;
         stEn(i,2) = nPnts(idx)-100;
         moveOn = 0;
      elseif boutMask(stEn(i,2))==1
         pIdx = find(pInf.wc == stEn(i,2));
         if pInf.ampR(pIdx)<6 && pIdx > 1 && ...
               (pInf.wc(pIdx-1)>stEn(i,1)) && pInf.ampR(pIdx-1)>20
            pul2rem(pIdx) = 1;
            moveOn = 0;
         end
      end
      if boutMask(stEn(i,1))==1
         pIdx = find(pInf.wc == stEn(i,1));
         if pInf.ampR(pIdx)<6 && pIdx < length(pInf.wc) && ...
               (pInf.wc(pIdx+1)<stEn(i,2) && pInf.ampR(pIdx+1)>20)
            pul2rem(pIdx) = 1;
            moveOn = 0;
         end
      end
   end
   [pInf, pul2rem] = remBadPul(pInf, pul2rem);
   [stEn, boutMask] = findSongBouts(pInf, wInf, nAmp, blkP, 1);
end
[wInf] = combineSineSong(sInf, pInf, wInf, blkP, 1);
[bInf.stEn, bInf.Mask] = findSongBouts(pInf, wInf, nAmp, blkP, 1);
pul2rem(bInf.Mask(pInf.wc) == 0) = 1;
[pInf] = remBadPul(pInf, pul2rem);
[~, pIdx] = setdiff(pInfOrig.wc, pInf.wc);
for i = 1:length(pIdx)
   newWC = pInfOrig.wc(pIdx(i));
   if bInf.Mask(newWC) > 0 && wInf.aSin(newWC) == 0 &&...
         min(abs(pInf.wc-newWC)) > 200 && pInfOrig.ampR(pIdx(i)) > 3;
      idx = length(pInf.wc)+1;
      pInf.wc(idx) = pInfOrig.wc(pIdx(i));
      pInf.detM(idx,:) = pInfOrig.detM(pIdx(i),:);
      pInf.mAmp.Mx(idx,:) = pInfOrig.mAmp.Mx(pIdx(i),:);
      pInf.mAmp.Int(idx,:) = pInfOrig.mAmp.Int(pIdx(i),:);
      pInf.pSec(idx,:) = pInfOrig.pSec(pIdx(i),:);
      pInf.pLik(idx,:) = pInfOrig.pLik(pIdx(i),:);
      pInf.hpWi(idx,:) = pInfOrig.hpWi(pIdx(i),:);
      pInf.c2SR(idx,:) =  pInfOrig.c2SR(pIdx(i),:);
      pInf.ampR(idx,:) =  pInfOrig.ampR(pIdx(i),:);
   end
end
pul2rem = 0*pInf.wc';
pInf = remBadPul(pInf, pul2rem);
wInf = combineSineSong(sInf,pInf,wInf,blkP,bInf.Mask);
bInf = finalizeBoutInfo(pInf, wInf, bInf);

% finalize everything
Song = int16(pInf.Song*dataScalingFactor);
pInf.pSec = int16(pInf.pSec*dataScalingFactor);
pInf = rmfield(pInf, {'Song'; 'ampR'; 'exSong'});
wInf = rmfield(wInf, {'nAmp'; 'bSin'});
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [detM, wc, c2SR] = findUniPul(sInf, oneSong)
sampLength = length(oneSong);
for i = 1:length(sInf)
   if isfield(sInf(i).pulseInfo2, 'scmx')
      sInf(i).pulseInfo2 = rmfield(sInf(i).pulseInfo2, 'scmx');
   end
end
for i = 1:length(sInf)
   if ~isfield(sInf(i).pulseInfo2, 'wc'); continue; end
   tempWC{i,1} = sInf(i).pulseInfo2.wc(sInf(i).pulseInfo2.wc>0);
   if isfield(sInf(i).pulseInfo2, 'c2SR')
      tempc2SR{i,1} = sInf(i).pulseInfo2.c2SR(sInf(i).pulseInfo2.wc>0);
   else
      tempc2SR{i,1} = [];
   end
end
pWid = 251;
aPul = [];
for i = 1:length(sInf)
   aPul = [aPul;[i*ones(length(tempWC{i}),1), tempWC{i}', tempc2SR{i}']];
end
[wcSrt, mIdx] = sort(aPul(:,2));
mIdx(wcSrt<ceil(pWid/2)) = [];
wcSrt(wcSrt<ceil(pWid/2)) = [];
mIdx(wcSrt>(sampLength-ceil(pWid/2))) = [];
wcSrt(wcSrt>(sampLength-ceil(pWid/2))) = [];
mIdx = aPul(mIdx,[1 3]);
dPnts = diff([0;wcSrt]);
sPnts = find(dPnts>50);
sPntConflict = find(dPnts>50 & dPnts<125,2);
while ~isempty(sPntConflict)
   probPnt = sPntConflict(1);
   sPntIdx = find(sPnts == probPnt);
   numPrePul = sPnts(sPntIdx) - sPnts(sPntIdx-1);
   if sPntIdx < length(sPnts)
      numPostPul = sPnts(sPntIdx+1) - sPnts(sPntIdx);
   else
      numPostPul = 100000;
   end
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
      wcSrt(probPnt:probPnt+numPostPul-1) = [];
      mIdx(probPnt:probPnt+numPostPul-1,:) = [];
   else
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
for i = 1:length(sPnts)
   if i == length(sPnts)
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
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wInf] = combineSineSong(sInf, pInf, wInf, blkP, bMsk)
oneSong = pInf.Song;
wc = pInf.wc;
nAmp = wInf.nAmp;
nThresh = 0.15;
if ~isfield(wInf, 'aSin')
   aSin = 0*oneSong;
   bSin = 0*oneSong;
   for i = 1:length(sInf)
      if ~isfield(sInf(i).winSine, 'NWK3060') || isempty(sInf(i).winSine.NWK3060)
         sInf(i).winSine.NWK3060.num_events = [];
         sInf(i).winSine.NWK3060.start = [];
         sInf(i).winSine.NWK3060.stop = [];
         sInf(i).winSine.NWK3060.length = [];
      end
      if ~isfield(sInf(i).winSine, 'NWK2012') || isempty(sInf(i).winSine.NWK2012)
         sInf(i).winSine.NWK2012.num_events = [];
         sInf(i).winSine.NWK2012.start = [];
         sInf(i).winSine.NWK2012.stop = [];
         sInf(i).winSine.NWK2012.length = [];
      end
      
      for j = 1:size(sInf(i).winSine.NWK3060.start,1)
         sRng = round(sInf(i).winSine.NWK3060.start(j):...
            sInf(i).winSine.NWK3060.stop(j));
         aSin(sRng) = aSin(sRng)+1;
      end
      
      for j = 1:size(sInf(i).winSine.NWK2012.start,1)
         sRng = round(sInf(i).winSine.NWK2012.start(j):...
            sInf(i).winSine.NWK2012.stop(j));
         bSin(sRng) = bSin(sRng)+1;
      end
   end
   bTest = smooth(double(abs(oneSong)>nAmp).*(bSin>0),1001)>0.15;
   wInf.aSin = 2*aSin + bTest';
   wInf.bSin = bSin;
   tDat=[false, wInf.aSin~=0, false];
   stEn = [strfind(tDat,[0,1]); strfind(tDat,[1,0])-1]';
end
if exist('blkP', 'var') && blkP(1)==5
   for i = 1:size(wInf.stEn,1)
      if sum(abs(oneSong(wInf.stEn(i,1):wInf.stEn(i,2))) > nThresh) > 50
         wInf.aSin(wInf.stEn(i,1):wInf.stEn(i,2)) = ...
            wInf.bSin(wInf.stEn(i,1):wInf.stEn(i,2));
      end
   end
end
tDat=[false, wInf.aSin~=0, false];
stEn = [strfind(tDat,[0,1]); strfind(tDat,[1,0])-1]';
i = 1;
while i ~= size(stEn,1) && size(stEn,1)~=0 && exist('blkP', 'var') && ...
      blkP(1) == 1;
   nBPnt = blkP(find(blkP>stEn(i,2),1));
   if any(wc<stEn(i+1,1)&wc>stEn(i,2)) || nBPnt<stEn(i+1,1); i = i+1;
      continue; end
   if stEn(i+1,1)-stEn(i,2)<1500
      wInf.aSin(stEn(i,2):stEn(i+1,1)) = 1;
      stEn(i,2) = stEn(i+1,2);
      stEn(i+1,:) = [];
   else
      i = i+1;
   end
end
dIdx = [];
for i = 1:size(stEn,1)
   sRng = stEn(i,1):stEn(i,2);
   mTest = mean(abs(oneSong(sRng)))./std(abs(oneSong(sRng)));
   sTest = mean(abs(oneSong(sRng))>nAmp);
   oneSong(sRng) = oneSong(sRng).*...
      abs(oneSong(sRng)) < 5*mean(abs(oneSong(sRng)));
   numDet = mean(wInf.aSin(sRng));
   if diff(stEn(i,:)) < 700 || (diff(stEn(i,:))<1000 && numDet<2.5)...
         || (sTest<0.15 && numDet<3) || numDet == 1 || ...
         (mTest < 1 && numDet<3 && exist('bMsk', 'var'));
      dIdx = [i;dIdx];
      wInf.aSin(stEn(i,1):stEn(i,2)) = 0;
   end
end
stEn(dIdx,:) = [];
if exist('bMsk', 'var') &&  sum(bMsk(:))~=1
   wInf.aSin(bMsk == 0) = 0;
   for i = 1:length(wc)
      wInf.aSin(wc(i)+[-50:50]) = 0;
   end
   tDat=[false,wInf.aSin~=0,false];
   stEn = [strfind(tDat,[0,1]); strfind(tDat,[1,0])-1]';
end
wInf.stEn = stEn;
wInf.aSin = cast(wInf.aSin, 'int8');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf,pul2rem] = remBadPul(pInf, pul2rem)
pInf.detM(pul2rem>0,:) = [];
pInf.wc(pul2rem>0,:) = [];
pInf.mAmp.Mx(pul2rem>0,:) = [];
pInf.mAmp.Int(pul2rem>0,:) = [];
pInf.pSec(pul2rem>0,:) = [];
pInf.hpWi(pul2rem>0,:) = [];
pInf.c2SR(pul2rem>0,:) = [];
pInf.ampR(pul2rem>0,:) = [];
if isfield(pInf, 'pLik')
   pInf.pLik(pul2rem>0) = [];
end
pul2rem = 0*pInf.wc;

[~, srtIdx] = sort(pInf.wc);
pInf.detM = pInf.detM(srtIdx,:);
pInf.wc = pInf.wc(srtIdx,:);
pInf.mAmp.Mx = pInf.mAmp.Mx(srtIdx,:);
pInf.mAmp.Int = pInf.mAmp.Int(srtIdx,:);
pInf.pSec = pInf.pSec(srtIdx,:);
pInf.hpWi = pInf.hpWi(srtIdx,:);
pInf.c2SR = pInf.c2SR(srtIdx,:);
pInf.ampR = pInf.ampR(srtIdx,:);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stEn, boutMask] = findSongBouts(pInf, wInf, nAmp, blkP, fin)
sBouts = pInf.Song*0;
sBouts(pInf.wc) = 1;
sBouts(wInf.aSin>0) = 1;
tDat = [false,sBouts~=0,false];
stEn = [strfind(tDat,[0,1]); strfind(tDat,[1,0])-1]';
sBouts(pInf.wc) = 10;
sBouts(wInf.aSin>0) = 2;
idx = 1;
while idx < size(stEn,1)
   nBPnt = blkP(find(blkP>stEn(idx,2),1));
   if stEn(idx+1,1)-stEn(idx,2)<2000 && nBPnt>stEn(idx+1,1)
      sBouts(stEn(idx,2)+1:stEn(idx+1,1)-1) = 1;
      stEn(idx,:) = [stEn(idx,1) stEn(idx+1,2)];
      stEn(idx+1,:) = [];
   else
      idx = idx + 1;
   end
end
stEn(diff(stEn, [], 2)<700,:) = [];  % merge IBIs<70ms
% assemble bouts to delete
dIdx = [];  % will contain idx of bouts to delete
for i = 1%:size(stEn,1)
   bRng = stEn(i,1):stEn(i,2);
   nEvs = [sum(sBouts(bRng)==10), (sum(sBouts(bRng)==2))/1000]; % count number pulses (`==10`) and milliseconds of sine song (`==2`)
   pIdx = pInf.wc>=bRng(1)&pInf.wc<=bRng(end);     % idx to get song trace
   snRng = pInf.Song(bRng);                        % song trace for that bout
   % FIX: this deletes the first 3.5 seconds - we don't want that!!!
   %    if bRng(1)<35000 || bRng(end)>length(pInf.Song)-35000 || ...      % if in the first or last 3.5seconds
   if (nEvs(1)<3 && nEvs(2)<1) || (nEvs(1)<2 && nEvs(2)<2)        % or <3 pulses and <1ms sine or <2 pulses and <2ms sine
      dIdx = [dIdx;i];
   elseif ((mean(pInf.pLik(pIdx)>10)<=0.5 || mean(pInf.pLik(pIdx)>0)<0.6) ...  % if unlikely pulses
         && (nEvs(1)<6 && nEvs(2)<2)) && mean(pInf.ampR(pIdx)>3)>=0.4          % and <6pulses and <2ms sine and soft pulses
      dIdx = [dIdx;i];
   elseif nEvs(2) == 0 && exist('fin', 'var')                           % if no sine - pulse only bout
      pTest = diff(find(sBouts(bRng)==10)); % ipis for that bout
      if (sum(pTest>250 & pTest<750)<3 && ...      %% if <3 "normal" ipis
            (any(pInf.pLik(pIdx)<10) || any(pInf.hpWi(pIdx)<0.6)))...   %  and unlikely or short pulses and
            || (nEvs(1)<4 && mean(pInf.ampR(pIdx)>5)<1)                 % or if few and soft pulses
         dIdx = [dIdx;i];
      end
   elseif exist('fin', 'var') && ((nEvs(1)<3 && nEvs(2)<1) || ... % if <3pulses and <1ms sine
         (nEvs(1) == 0 && nEvs(2)<3))                             % or if no pulses and <3ms sine
      dIdx = [dIdx;i];
   elseif nEvs(1)<3 && exist('fin', 'var')
      sTest = mean(abs(snRng(sBouts(bRng) == 2))>nAmp);  % sine song for that bout that is above noise floor
      if (sTest<0.2 && (mean(wInf.aSin(bRng))<3 || ...   % if too little sine song above noise and little sinse song(?)
            (nEvs(2)<5 && nEvs(1)<2))) || sTest<0.1
         dIdx = [dIdx;i];
      end
   end
end
stEn(dIdx,:) = [];
boutMask = int8(pInf.Song*0 - 1);
for i = 1:size(stEn,1)
   boutMask(stEn(i,1):stEn(i,2)) = 1;
end
boutMask(wInf.aSin>0) = boutMask(wInf.aSin>0)+1;
boutMask(boutMask==-1) = 0;
if exist('fin', 'var')
   idx = 1;
   while idx < size(stEn,1)
      nBPnt = blkP(find(blkP>stEn(idx,2),1));
      if stEn(idx+1,1)-stEn(idx,2) < 2000 && ...
            nBPnt>stEn(idx+1,1)
         boutMask(stEn(idx,2)+1:stEn(idx+1,1)-1) = 1;
         stEn(idx,:) = [stEn(idx,1) stEn(idx+1,2)];
         stEn(idx+1,:) = [];
      else
         idx = idx + 1;
      end
   end
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bInf] = finalizeBoutInfo(pInf, wInf, bInf)
if isempty(bInf.stEn)
   bInf.stEn = [];
   bInf.Mask = [];
   bInf.Buff = [];
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

buff = [35000 35000]; % why is this hard-coded?
tempSin = find(wInf.aSin)';
temp = 0*bInf.Mask;
temp(pInf.wc) = 1;
temp(tempSin) = 2;
iPnts =  find(temp>0);
tDat = interp1(iPnts, double(temp(iPnts)),...
   1:length(temp), 'nearest');
bInf.Mask = tDat.*double(bInf.Mask>0);
bInf.Mask = int8(bInf.Mask);
Song = pInf.Song;
bInf.Buff = buff;
for i = 1:size(bInf.stEn,1)
   bInf.Bout{i,1} = Song(max(bInf.stEn(i,1)-buff(1),1):min(length(Song), bInf.stEn(i,2)+buff(2)));
   % pre and postpend nan if buffered bout exceeds song duration (<0, >length(Song)
   bInf.Bout{i,1} = [nan(1, -min(bInf.stEn(i,1)-buff(1),0)), bInf.Bout{i,1}, nan(1, -min(length(Song) - (bInf.stEn(i,2)+buff(2)),0))];
   if max(bInf.Bout{i,1})<3
      bInf.Bout{i,1} = int16(bInf.Bout{i,1}*1000);
   end
   bInf.bRng(i,:) = [buff(1)+1 (length(bInf.Bout{i}) - buff(2))];
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
function [pInf] = findSongNearPulse(pInf,wInf)
pInf.proxSng = 0*pInf.wc;
tempSin = find(wInf.aSin)';
if numel(pInf.wc) > 1
   for i = 1:numel(pInf.wc)
      if i == 1
         tstPul = [pInf.wc(i+1); tempSin(1:200:end)];
      elseif i == numel(pInf.wc)
         tstPul = [pInf.wc(i-1); tempSin(1:200:end)];
      else
         tstPul = [pInf.wc(i-1); pInf.wc(i+1); tempSin(1:200:end)];
      end
      tPul = abs(tstPul - pInf.wc(i));
      pInf.proxSng(i,1) = min(tPul);
   end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf, pul2rem] = remWorstPul(pInf,wInf,pul2rem)
moveOn = 0;
while moveOn == 0
   [pInf,pul2rem] = remBadPul(pInf, pul2rem);
   pInf = findSongNearPulse(pInf,wInf);
   moveOn = 1;
   for i = 1:length(pInf.wc)
      [nDis, nIdx] = sort(abs(pInf.wc - pInf.wc(i)));
      amp = [pInf.ampR(i); pInf.ampR(nIdx(2:3))];
      type = length(setdiff([i-1 i+1], nIdx(2:3)));
      type = 1+(type - 1)*-1;
      % no idea what this is doing!!
      if nDis(2)<150 && amp(2)>amp(1) || ...
            (type==1 && amp(2)>5*amp(1) && nDis(2)<200) || ...
            (type==1 && amp(2)>7*amp(1) && amp(1)<2.5) || ...
            (type==2 && min(amp(2:3))>3*amp(1) && min(nDis(2:3))<300) || ...
            (type==2 && amp(1)<1.75 && min(nDis(2:3))<250) || ...
            (sum(pInf.detM(i,:))==1 && amp(1)<1.75) || ...
            (sum(nDis(2:3)<800)<type && (amp(1)<2.5 || pInf.hpWi(i)<0.5))|| ...
            (amp(1)<10 && pInf.proxSng(i)>1200) || pInf.hpWi(i)<0.3 ...
            || pInf.c2SR(i) < 0.25;
         pul2rem(i) = 1;
         moveOn = 0;
      end
   end
end
[pInf,pul2rem] = remBadPul(pInf, pul2rem);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pInf, wInf, pul2rem] = remSinConPul(pInf, wInf, sInf, pul2rem)
sinPul = find(wInf.aSin(pInf.wc) > 0);
moveOn = 0;
nAmp = wInf.nAmp;
while moveOn == 0
   moveOn = 1;
   for i = 1:length(sinPul)
      tChk = pul2rem==0;
      if pInf.c2SR(sinPul(i))<2.5 && pul2rem(sinPul(i)) == 0;
         [dPul] = sort(abs(pInf.wc(sinPul(i))-pInf.wc(tChk)));
         if dPul(2) > 300;
            pul2rem(sinPul(i)) = 1;
            moveOn = 0;
         else
            sidAmps =[max(abs(pInf.pSec(sinPul(i),1:50))) ...
               max(abs(pInf.pSec(sinPul(i),201:251)))];
            pInf.c2SR(sinPul(i),1) = max(abs(pInf.pSec(sinPul(i),90:160)))...
               /max(sidAmps);
            if  pInf.c2SR(sinPul(i)) < 2.5
               pul2rem(sinPul(i)) = 1;
               moveOn = 0;
            end
         end
      end
      if (pInf.c2SR(sinPul(i))<4.5 || pInf.hpWi(sinPul(i)) < 0.7)...
            && pul2rem(sinPul(i)) == 0
         [dPul] = sort(abs(pInf.wc(sinPul(i))-pInf.wc(tChk)));
         if dPul(2)> 800 || (pInf.hpWi(sinPul(i)) < 0.6 && ...
               pInf.ampR(sinPul(i)) > 20)
            pul2rem(sinPul(i)) = 1;
            moveOn = 0;
         end
      end
   end
end
for i = 1:length(sinPul)
   if pul2rem(sinPul(i)) == 0
      wInf.aSin(pInf.wc(sinPul(i))-125:pInf.wc(sinPul(i))+125) = 0;
   end
end
[wInf] = combineSineSong(sInf, pInf,wInf);
moveOn = 0;
while moveOn == 0 && ~isempty(wInf.stEn)
   moveOn = 1;
   tIdx = find(pul2rem);
   for i = 1:length(tIdx)
      tstPul = pInf.wc(tIdx(i));
      if (min(abs(tstPul-wInf.stEn(:))) < 800 && pInf.c2SR(tIdx(i))>1.25) ...
            || wInf.aSin(tstPul) == 0
         sIdx = find(wInf.stEn(:,2)>=pInf.wc(tIdx(i)),1);
         [dPul] = sort(abs(pInf.wc(pul2rem==0)-tstPul));
         if wInf.aSin(tstPul) == 0 || dPul(1) < 600 || ...
               (dPul(1)>600 && dPul(1)<800 && (pInf.ampR(tIdx(i))>6 ...
               || diff(wInf.stEn(sIdx,:) < 1000)))
            if wInf.aSin(tstPul) == 0 || diff(wInf.stEn(sIdx,:))<1500 ...
                  || mean(pInf.exSong(wInf.stEn(sIdx,1):wInf. ...
                  stEn(sIdx,2)))/nAmp < 0.5*pInf.ampR(tIdx(i))
               pul2rem(tIdx(i)) = 0;
               wInf.aSin(tstPul-125:tstPul+125) = 0;
               moveOn = 0;
            end
         end
      end
   end
   [wInf] = combineSineSong(sInf, pInf,wInf);
end
[wInf] = combineSineSong(sInf, pInf,wInf, 5);
[pInf,pul2rem] = remBadPul(pInf, pul2rem);
end



