function FlatMasker( Signal_Length, SYS_or_LUT_resolution, Noise_Mask_dB, weight, setup )
%FLATMASKER Generates loudspeaker signals for a spatially flat noise masker

SYS_type = 'Current_Systems.SR_System';

%% Setup Variables
if nargin == 2
    if isa(SYS_or_LUT_resolution,SYS_type)
        SYS = SYS_or_LUT_resolution;        
        signal_info = SYS.signal_info;
        system_info = SYS.system_info;
        setup = SYS.Masker_Setup;
        if any(size(signal_info.L_noise_mask)~=1)
           warning('Different number of noise levels may not be supported.'); 
        end
    else
        error(['Second input argument must be of type: ' SYS_type]);
    end
elseif nargin == 5
    system_info.Drive = 'Z:\'; % Database drive (storage drive)
    system_info.LUT_resolution = SYS_or_LUT_resolution;
    
    signal_info.c = 343; % Speed of sound in metres/sec
    signal_info.Fs = 16000; % Sampling frequency
    signal_info.Nfft = 1024;% Number of fft components
    signal_info.overlap = 0.5;
    signal_info.f_low  = 150;  % Hz
    signal_info.f_high = 8000; % Hz
    signal_info.L_noise_mask = Noise_Mask_dB; % dB
    signal_info.weight = weight;
    signal_info.method = 'FlatMasker';
end

signal_info.input_filename = 'Masker';

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, system_info.LUT_resolution, system_info.Drive, 'new');



%% Generate Loudspeaker Signals
Loudspeaker_Signals = reshape( v_addnoise( ... % White noise
                            zeros(Signal_Length,setup.Loudspeaker_Count), ...
                            signal_info.Fs, -Inf), ...
                            [],setup.Loudspeaker_Count) ...
                            * db2mag(signal_info.clipLevelAdjust) ... % Adjust level to avoiding clipping
                            * db2mag(signal_info.L_noise_mask) ... % Adjust level to given value
                            / sqrt(setup.Loudspeaker_Count) ... % Adjust level for equal power
                            * 2 ; %<- Find out why this helps

                        
%% Once we have the speaker signals we should save them for later use as .wav files
Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    [], [], [], [], ...
    signal_info.Fs );


end


