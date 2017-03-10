function [nLevel, winSine, pulseInfo, pulseInfo2, pcndInfo, cmhSong] = segmentSong(xsong, varargin)
%USAGE
%  [noiseLevel, winnowed_sine, pulseInfo, pulseInfo2, pcndInfo] = Process_SongRondo(xsong, [params_path])
%PARAMS
%  xsong - single channel song recording
%  params_path - OPTIONAL path to parameter file that will get executed
%RETURNS:
%  [noiseLevel, winnowed_sine, pulseInfo, pulseInfo2, pcndInfo]
%

if nargin==2
   params_path = varargin{1};
else
   params_path = [];
end
FetchParams% set default parameters

disp('Running multitaper analysis on signal.')
disp(['Song duration is ' num2str(length(xsong)/Params.Fs/60,3) ' minutes.']);

if Params.find_sine
   disp('finding sine.')
   % first pass
   Sines.MultiTaper = MultiTaperFTest(xsong, Params.Fs, Params.NW, ...
      Params.K, Params.dT, Params.dS, Params.pval, Params.fwindow);
   winSine.NWK3060 = ...
      SineSegmenter(xsong, Sines.MultiTaper, Params.Fs, Params.dT, Params.dS, ...
      Params.sine_low_freq, Params.sine_high_freq, Params.sine_range_percent);
   
   % second pass - with Params.dT=12 Params.dS=23
   Sines.MultiTaper = MultiTaperFTest(xsong, Params.Fs, 12, ...
      23, Params.dT, Params.dS, Params.pval, Params.fwindow);
   winSine.NWK2012 = ...
      SineSegmenter(xsong, Sines.MultiTaper, Params.Fs, Params.dT, Params.dS, ...
      Params.sine_low_freq, Params.sine_high_freq, Params.sine_range_percent);
else
   disp('you chose not to find sine.')
   winSine.NWK3060 = [];
   winSine.NWK2012 = [];
   %%
   disp('doing it anyways for noise estimation')
   % we actually only need:
   % Sines.MultiTaper.f and Sines.MultiTaper.A for EstimateNoise
   % but this is relatively fast so doesn't hurt to run it
   Sines.MultiTaper = MultiTaperFTest(xsong, Params.Fs, 12, ...
      23, Params.dT, Params.dS, Params.pval, Params.fwindow);
end

fprintf('Finding noise floor.\n')
noise = EstimateNoise(xsong, Sines.MultiTaper, Params, Params.low_freq_cutoff,...
   Params.high_freq_cutoff);

disp('finding pulse.')
[pcndInfo, pulseInfo, pulseInfo2, cmhSong] = ...
   SegmentPulse(xsong,noise.d,Params.fc,Params.DoGwvlt,...
   Params.pWid,Params.minIPI,Params.thresh,Params.minAmplitude,...
   Params.maxIPI,Params.frequency,Params.close,Params.Fs);
nLevel = mean(abs(noise.d));

