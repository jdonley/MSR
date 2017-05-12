function [filts, IRs] = getCalibrationFilters( recording, referenceSig, SYS )

F_band = [SYS.system_info.f_low, SYS.system_info.f_high];
Fs = SYS.system_info.fs;
FiltLen = SYS.system_info.Calibration_FiltLen;
Reg = SYS.system_info.Calibration_FiltReg;

isLineArray = contains(SYS.system_info.CurrentSpeakerArrayType, 'line');
if isLineArray
    TFAlign = Speaker_Setup.Calibration.linTFalign( SYS );
end

Dur = size(referenceSig,1) / Fs; %seconds
Nspkrs = floor(size(recording,1) / Dur / Fs );

IRs = []; invIRs = [];
for chan = flip( 1:size(recording,2) ) % the flip allows 3D matrix initialisation by starting at the larger index
    
    % (commented out as using more memory friendly for loops instead)
    %Segment loudspeaker recordings
    %     [R, referenceSigNoPad] = Speaker_Setup.Calibration.SplitRecording( ...
    %         recording(:,chan), referenceSig, Nspkrs, Dur, Fs );
    
    %Determine impulse responses
    IRs(:,:,chan) = Speaker_Setup.Calibration.IRsFromRecordings( ...
        squeeze(recording(:,chan,:)), referenceSig, F_band, Fs );
    
    %Compute inverse impulse responses
    invIRs(:,:,chan) = Speaker_Setup.Calibration.getInverseFilters( ...
        IRs(:,:,chan), FiltLen, F_band, Fs, Reg, 'kirkeby' );
    
    if isLineArray
        invIRs(:,:,chan) = Tools.fconv( invIRs(1:end-size(TFAlign,2)+1,:,chan), TFAlign.' );
    end
end

filts = invIRs;

end