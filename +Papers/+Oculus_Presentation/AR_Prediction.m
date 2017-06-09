
SpeechFilePath = '\+Miscellaneous\+Speech_Files\Female_SA2.WAV';


sig_info.Nfft = 12 * 1e-3 * sig_info.Fs; % Number of fft components
sig_info.time_delay = 0; % Seconds %If empty the time delay will based on the frame length
sig_info.overlap = 0.50;
sig_info.zeropadtime = sig_info.Nfft / sig_info.Fs; % miliseconds of zero padding (i.e. maximum time shift before circular convolution)
sig_info.predict_buff = 2; % Length of prediction buffer as a percentage of a full frame prediction (i.e. 10 is %1000 of frame length)
sig_info.AR_method = 'lpc';
sig_info.window = false;



[ ~, ~, Frames, ~, ~, Frames_orig] = ...
    Broadband_Tools.FFT_custom(SpeechFilePath, sig_info);
