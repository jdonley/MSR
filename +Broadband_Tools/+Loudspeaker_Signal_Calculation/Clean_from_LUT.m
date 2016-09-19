function Clean_from_LUT( Input_file, SYS_or_LUT_resolution, weight, setup )
%CLEAN_FROM_LUT Generates spatially weighted loudspeaker signals for the
%reproduction of zoned clean speech.

SYS_type = 'Current_Systems.SR_System';

%% Setup Variables
if nargin == 2
    if isa(SYS_or_LUT_resolution,SYS_type)
        SYS = SYS_or_LUT_resolution;        
        signal_info = SYS.signal_info;
        system_info = SYS.system_info;
        setup = SYS.Main_Setup;
    else
        error(['Second input argument must be of type: ' SYS_type]);
    end
elseif nargin == 4
    system_info.Drive = 'Z:\';
    system_info.LUT_resolution = SYS_or_LUT_resolution;
    
    signal_info.c = 343; % Speed of sound in metres/sec
    signal_info.Fs = 16000; % Sampling frequency
    signal_info.Nfft = 1024;% Number of fft components
    signal_info.overlap = 0.5;
    signal_info.f_low  = 150;  % Hz
    signal_info.f_high = 8000; % Hz
    signal_info.f_high_meas = 7000; % Hz %Maximum frequency with accurate response at given sampling rate
    signal_info.weight = weight;
    signal_info.method = 'NoMask';
end

% Clean speech means no noise levels
signal_info.L_noise_mask = -Inf; % dB
    
[~, Input_file_name, Input_file_ext] = fileparts( Input_file );
signal_info.input_filename = Input_file_name;



[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, system_info.LUT_resolution, system_info.Drive, 'new');



%% Read input signal
Input_Signal = audioread( Input_file );

if ~isempty(strfind(lower(signal_info.method), 'cancel'))
    %%% TODO
    % Fix this magnitude adjustment which is dependent on the transfer
    % function for point sources. Originates somewhere in
    % loudspeaker_setup.m or earlier.
    magnitudeADJ = 0.8;
    Input_Signal = conv(-1,Input_Signal) * magnitudeADJ;
end

%% Generate Loudspeaker Signals
[Loudspeaker_Signals, Original_] = Broadband_Tools.Loudspeaker_Signal_Calculation.getMSRLoudspeakerSignals( ...
    Input_Signal, setup, signal_info, system_info );

%% Normalise Loudspeaker Signals
% scaler = rms(Loudspeaker_Signals(:));
% Loudspeaker_Signals = Loudspeaker_Signals ./ scaler .* db2mag(-35)  ; % Hope -35dB RMS level is enough to avoid clipping when saving
Original_ = Original_ ./ max(abs(Original_(:)));




%% Once we have the speaker signals we should save them for later use as .wav files

Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    setup.Loudspeaker_Count, ...
    Original_, ...
    signal_info.input_filename, ...
    Input_file_ext, ...
    signal_info.Fs );



end