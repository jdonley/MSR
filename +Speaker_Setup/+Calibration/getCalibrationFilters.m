function [filts, IRs] = getCalibrationFilters( recording, referenceSig, F_band, Fs, FiltLen, Reg)

Dur = size(referenceSig,1) / Fs; %seconds
Nspkrs = floor(size(recording,1) / Dur / Fs );

IRs = []; invIRs = [];
for chan = flip( 1:size(recording,2) )
    
    %Segment loudspeaker recordings
%     [R, referenceSigNoPad] = Speaker_Setup.Calibration.SplitRecording( ...
%         recording(:,chan), referenceSig, Nspkrs, Dur, Fs );
    
    %Determine impulse responses
    IRs(:,:,chan) = Speaker_Setup.Calibration.IRsFromRecordings( ...
        squeeze(recording(:,chan,:)), referenceSig, F_band, Fs );
    
    %Compute inverse impulse responses
    invIRs(:,:,chan) = Speaker_Setup.Calibration.getInverseFilters( ...
        IRs(:,:,chan), FiltLen, F_band, Fs, Reg, 'kirkeby' );
    filts = invIRs;
    
end

end