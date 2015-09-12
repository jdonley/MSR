function [ Y, frequencies, Frames, Windows] = FFT_custom_vec( x, Nfft, Fs, Overlap )
%FFT_CUSTOM Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
   Nfft = 1024; % 1024 samples per frame
   Fs = 16000; % 16kHz sampling frequency
   Overlap = 0.5;
end

%% Read audio file
%[x, Fs_file] = audioread( audiofilepath );
if size(x,2) > 1 % If audio is not mono then sum audio from all channels together
   %x = sum(x,2);
end
x = x(:,1);

%% Split audio file into frames with overlap
% Find lengths of the overlapping and non-overlapping sections
non_overlap_size = floor(Nfft*(1-Overlap));
overlap_size = Nfft - non_overlap_size;

% Pad the signal file to allow the creation of the matrix
N_of_frames = ceil( length(x) / non_overlap_size );
r = rem( length(x), non_overlap_size);
r(r==0)=non_overlap_size;
pad = zeros( non_overlap_size - r , 1 );
x = [x; pad]';

% Find the overlapping section and non-overlapping section
first_part_frame  = reshape(x, [non_overlap_size N_of_frames])';
second_part_frame = [first_part_frame(2:end, 1:overlap_size); zeros(1,overlap_size)];

% Concatenate the two sections
Frames = [first_part_frame second_part_frame ];

% Remove last frame if it contains no extra information other than an overlap
if nnz(Frames(end,:)) <= overlap_size
    Frames = Frames(1:(end-1),:);
    N_of_frames = size(Frames,1);
end


%% Apply windowing to each frame
% if mod(Nfft,2) % If odd correct matlabs silly bugger default hamming window    
%     window_ = hamming(Nfft);
%     window_(1) = window_(1)/2;  % repair constant-overlap-add for 50% overlap
%     window_(Nfft) = window_(Nfft)/2;
% else
%     window_ = hamming(Nfft+1);
%     window_ = window_(1:Nfft); % Delete last sample for periodic case
% end

if mod(Nfft,2) % If odd correct matlabs silly bugger default hamming window    
    window_ = hanning(Nfft);
else
    window_ = [hanning(Nfft-1); 0];
end

%Square root window
window_ = sqrt(window_);

Windows = repmat( window_', [N_of_frames 1]);
Frames = Frames .* Windows;

%% Compute FFT of each frame
Y = fft(Frames, [], 2) / Nfft;
frequencies = linspace(0, Fs/2, Nfft/2+1);
frequencies = frequencies(1:end-1); % We do not return the Nyquist Frequency because it folds to the equivalent of the negative DC frequency

%% Return Amplitude and Phase for every frequency for each frame.
Y = Y(:,1:Nfft/2); % One-sided spectrum from (assumed) real signal input


end

