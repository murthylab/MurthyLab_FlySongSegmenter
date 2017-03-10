function oneSong = mergeChannels(data, doFilter)
% oneSong = mergeChannels(data)
% merge all recorded channels to a single trace using 'max'
if ~exist('doFilter', 'var')
   doFilter = true;
end

Fs = 1e4;% Hz
[duration, NumberOfChannels] = size(data);
if doFilter
   % init BP filter
   Hd = fdesign.bandpass('N,F3dB1,F3dB2',10,80,460,Fs);
   Hd = design(Hd,'butter');
   % filter data
   for chan = 1:NumberOfChannels
      data(:,chan) = filtfilt(Hd.sosMatrix,Hd.ScaleValues, data(:,chan));
   end
end
% get max value in +/-5 ms neighbourhoods
dataMax = runningExtreme(abs(data), 101, 'max');
% merged song is given by the maximum across channels
% not sure this is the best way to combine the data...
[~, mIdx] = max(dataMax, [], 2);
clear dataMax
% convert to linear indices so we can assign directly
ind = sub2ind([duration, NumberOfChannels], 1:duration, mIdx');
oneSong = data(ind);