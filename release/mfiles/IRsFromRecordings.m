function LoudspeakerIRs = IRsFromRecordings( SplitRecordingMatrix, InputTestSignal, f_band, fs, ir_duration )
%IRSFROMRECORDING Find the impulse response for each recording
if nargin < 5
    Nfft = length(InputTestSignal);
else
    Nfft = ir_duration*fs*2;
end

%A = load('+Tools\invsweepfft.mat');
invsweepfft = Tools.invSweepFFT( InputTestSignal, f_band(1), f_band(2), fs, Nfft);

Nspkrs = size(SplitRecordingMatrix,2);
LoudspeakerIRs = zeros( Nfft/2, Nspkrs );

for spkr = 1:Nspkrs
    LoudspeakerIRs( :, spkr ) = Tools.extractIR( ...
        SplitRecordingMatrix( :, spkr ), ...
        invsweepfft);
end

end

