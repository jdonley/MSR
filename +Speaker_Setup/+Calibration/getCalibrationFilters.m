function [filts, IRs] = getCalibrationFilters( recording, referenceSig, F_band, Fs, FiltLen, Reg)

Dur = length(referenceSig) / Fs; %seconds
Nspkrs = floor(length(recording) / Dur / Fs );

%Segment loudspeaker recordings
[R, referenceSigNoPad] = Speaker_Setup.Calibration.SplitRecording( recording, referenceSig, Nspkrs, Dur, Fs );

%Determine impulse responses
IRs = Speaker_Setup.Calibration.IRsFromRecordings( R, referenceSigNoPad, F_band, Fs );

%Compute inverse impulse responses
invIRs = Speaker_Setup.Calibration.getInverseFilters( IRs, FiltLen, F_band, Fs, Reg, 'kirkeby' );
filts = invIRs;

% %Equalise frequency response
% irs = Tools.fconv( IRs, invIRs );
% 
% %Compute gain adjustment
% IRgains = Speaker_Setup.Calibration.getGainAdjustments( irs, F_band, Fs );
% 
% %Compute time-alignment
% TimeAlignIRs = Speaker_Setup.Calibration.getTimeAlignIRs( IRs );
% 
% %Create filters
% filts = Tools.fconv( ...
%     invIRs .* repmat( IRgains, [size(invIRs,1), 1]), ...
%     TimeAlignIRs );

end