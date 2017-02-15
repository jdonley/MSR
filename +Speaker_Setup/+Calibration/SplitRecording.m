function [ Responses, Input_NoPad ] = SplitRecording( Recording, InputSignal, Nspkrs, Duration, Fs, PaddingDuration)
%SPLITRECORDING Splits a single recording file into a matrix containing
%each speakers individual response recoridng

if nargin < 6
    PaddingDuration = [0, 0];
end

% Remove extra samples from end
Recording( int32(Nspkrs*Duration*Fs+1):end ) = [];

% Split
Responses = reshape( ...
    Recording, ...
    int32(Duration * Fs), ...
    Nspkrs);

% Remove padding samples
Responses( ((Duration-PaddingDuration(2))*Fs+1):end, : ) = [];
Responses(                1:(PaddingDuration(1)*Fs), : ) = [];
InputSignal( ((Duration-PaddingDuration(2))*Fs+1):end ) = [];
InputSignal(                1:(PaddingDuration(1)*Fs) ) = [];
Input_NoPad = InputSignal;
end

