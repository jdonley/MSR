function [ Y, frequencies, Frames, Windows, Yo] = FFT_custom( audio_pathORdata, signal_info ) %Nfft, Fs, Overlap, pad, delay, BuffLen )
%FFT_CUSTOM Summary of this function goes here
%   Detailed explanation goes here
Nfft = signal_info.Nfft;
Fs = signal_info.Fs;
Overlap = signal_info.overlap;
pad = signal_info.zeropadtime;
delay = signal_info.time_delay;
BuffLen = signal_info.predict_buff;

%% Read audio file
in_type=whos('audio_pathORdata');in_type=in_type.class;
if strcmp( in_type, 'char' )
    
    [x, Fs_file] = audioread( audio_pathORdata );
    
    if Fs ~= Fs_file && ~isempty(Fs)
        error(['Sampling frequency of the file (' ...
            num2str(Fs_file) ...
            'Hz) does not match the given sampling frequency (' ...
            num2str(Fs_file) 'Hz) .']);
    end
    Fs = Fs_file;
else
    
    x = audio_pathORdata;
    
    if isempty(Fs)
        Fs = 16000; % 16kHz sampling frequency
    end
end

% if nargin < 2
%     Nfft = 1024; % 1024 samples per frame
% end
% if nargin < 4
%     Overlap = 0.5;
% end
% if nargin < 5
%     pad = [];
% end
% if nargin < 6
%     delay = 0;
% end
% if nargin < 7
%     BuffLen = 1;
% end
if isempty(delay)
    delay = Nfft / Fs; % Assume delay is equivalent to frame length
    %(i.e. how long it takes to get the frame is the number of samples in that frame)
end

if ~isempty(pad)
    Npad = pad * Fs;
%     Nframe = Nfft + Npad;
    %    Nfft = Nfft / 2;
end


if size(x,2) > 1 % If audio is not mono then sum audio from all channels together
    %x = sum(x,2);
end
x = x(:,1);

%% Split audio file into frames with overlap
Frames = enframe( x, Nfft, (1-Overlap)*Nfft );
N_of_frames = size(Frames,1);

%% Predict ahead and frame
Frames_orig = Frames;
 if signal_info.predict_buff ~= 0
    Inc = (1-Overlap)*Nfft;
    
    Predict_Buffer = [zeros( Nfft/Inc * (signal_info.predict_buff-1), Nfft*BuffLen); ...
        enframe( [x((Nfft-delay*Fs+1):end); zeros(Nfft-delay*Fs,1)], ...
        Nfft*BuffLen, Inc )];

    if isfield(signal_info,'AR_method')
        Frames = Broadband_Tools.PredictiveFraming( ...
            Frames, Predict_Buffer, delay*Fs , signal_info.AR_method );
    else
        Frames = Broadband_Tools.PredictiveFraming( Frames, Predict_Buffer, delay*Fs , 'burg' );
    end
 end

%% Zero pad
if ~isempty(pad)
    
    Frames_orig = [zeros(N_of_frames,Npad/2), ...
        Frames_orig, ...
        zeros(N_of_frames,Npad/2)];
    
    N_of_frames = size(Frames,1);
    Frames = [zeros(N_of_frames,Npad/2), ...
        Frames, ...
        zeros(N_of_frames,Npad/2)];
end


%% Apply windowing to each frame
if mod(Nfft,2) % If odd correct matlabs silly bugger default hamming window
    window_ = hanning(Nfft);
else
    window_ = [hanning(Nfft-1);0];
end

%Square root window
window_ = sqrt(window_);

% Hanning window sums to unity on overlap add for 50% overlap
% So adjust for more overlaps
IncF = 1/(1-Overlap); % 1 / Increment
if floor(log2(IncF))==log2(IncF) %if increment is pow of 2 less than 50% overlap
    window_ = window_ / (IncF/2);
end

Windows = repmat( window_.', [N_of_frames 1]);

if ~isempty(pad)
    Windows = [zeros(N_of_frames,Npad/2), ...
        Windows, ...
        zeros(N_of_frames,Npad/2)];
end

Frames = Frames .* Windows;
Frames_orig = Frames_orig .* Windows(1:size(Frames_orig,1),:);
%% Compute FFT of each frame
Y  = fft(Frames, [], 2) / size(Frames,2);
Yo = fft(Frames_orig, [], 2) / size(Frames_orig,2);
frequencies = linspace(0, Fs/2, size(Frames,2)/2+1);
frequencies = frequencies(1:end-1); % We do not return the Nyquist Frequency because it folds to the equivalent of the negative DC frequency

%% Return Amplitude and Phase for every frequency for each frame.
Y = Y(:,1:size(Frames,2)/2); % One-sided spectrum from (assumed) real signal input
Yo = Yo(:,1:size(Frames_orig,2)/2); % One-sided spectrum from (assumed) real signal input


end

