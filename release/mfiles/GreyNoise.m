function [ y ] = GreyNoise( T, Fs, level)
%This function returns an equal loudness noise
% 
% Syntax:	Y = GREYNOISE( T, FS, LEVEL)
% 
% Inputs: 
% 	T - Length of the grey noise in seconds
% 	Fs - Sampling frequency of the grey noise output
% 	level - Magnitude gain/level adjustment as a multiplier
% 
% Outputs: 
% 	y - Grey noise output signal vector
% 
% Example: 
% 	[ y ] = GreyNoise( 10, 16000, 0.8); % 10 seconds of grey noise at 16kHz sampling frequency and 0.8 gain.
% 
% See also: awgn, randn, wgn

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015
% Date: 10 October 2015 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <3
    level = 1;
end

N = T * Fs;
if rem(N,2)
    M = N+1;
else
    M = N;
end


x = randn(1, M) / (2*pi); % WGN

% Frequency Domain Tranform
X = fft(x); % Fast Fourier Transform

% Frequencies
NumPts = M/2 + 1;
freqs = linspace(0, Fs/2, NumPts);
cutoff = 20; %Hz

% Equal Loudness Levels
EqLevel_dB = Perceptual_Tools.Threshold_in_Quiet(freqs(freqs>cutoff),'ISO226')';
EqLevel_dB = [linspace(0,EqLevel_dB(1), length(freqs)-length(EqLevel_dB)), EqLevel_dB];

% Apply magnitude weighting
X(1:NumPts) = X(1:NumPts) .* db2mag(EqLevel_dB);

% Apply conjugation for negative frequency side of spectrum
X(NumPts+1:M) = conj(X(M/2:-1:2));

% Time Domain Transform
y = ifft(X); % Inverse Fast Fourier Transform

% prepare output vector y
y = real(y(1, 1:N));

% ensure unity standard deviation and zero mean value
y = y - mean(y);
yrms = sqrt(mean(y.^2));
y = y/yrms;
% Level adjustment
y = y * level;

y = applyTaperWindow(y, Fs);
end

function y = applyTaperWindow(y, Fs)
rampLen = 0.1; %seconds
zerosLen = 0.05; %seconds

% % 150ms tapering window
window_ = [zeros(1,zerosLen*Fs), linspace(0,1,rampLen*Fs)];
winLen = length(window_);

% Apply to beinging and end of signal
y(1:winLen) = y(1:winLen) .* window_;
y((end-winLen+1):end) = y((end-winLen+1):end) .* flip(window_);
end
