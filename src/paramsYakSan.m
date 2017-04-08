% custom parameters for segmenting YAKUBA and SANTOMEA song

Params.maxIPI = round(Params.Fs/2);   % if no other pulse within this many ticks, do not count as a pulse
Params.frequency = 1500;              % if best matched scale is greater than this frequency, do not count as a pulse

Params.find_sine = false;              % yakuba and santomea do not produce sine song - makes song segmentation MUCH faster

%estimating noise - yak/san song has much higher frequency
Params.low_freq_cutoff = 200;         % exclude data below this frequency
Params.high_freq_cutoff = 800;        % exclude data above this frequency
Params.cutoff_sd = 2;                 % exclude data above this multiple of the std. deviation - more restrictive now

