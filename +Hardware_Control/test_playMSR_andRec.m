clc;
clear;

%%
array_type = 'circle';
spkr_radius = 1.3;
N_spkrs = 24;

system_info.dev_model = 'ASIO Hammerfall DSP';
system_info.fs = 48000;
system_info.playbackChannels = ...
    [ 1  2  3  4  5  6  7  8 ...
    9 10 11 12 13 14 15 16 ...
    17 18 19 20 21 22 23 24];
% First half recordings in Bright zone, second half in Quiet zone
system_info.recordChannels = ...
    [ 1 2 ];

system_info.sc = '_'; % Separating character for ascii paths
system_info.Drive = 'Z:\';
system_info.Filter_dir = '+Speaker_Signals\CalibratedSweeps\';
system_info.Calibrated_Signals_dir = '+Calibrated_Speaker_Signals\';

system_info.LUT_resolution = '512f_32w';

signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency of loudspeaker signals
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.L_noise_mask = -Inf; % dB
signal_info.weight = 1e4; %this is actually auto-calculated from maximum contrast
signal_info.method = 'NoMask';
signal_info.input_filename = [];

if strcmp(array_type, 'circle')
    [x,y] = pol2cart(-90/180*pi, 0.6);
    x_ = sqrt(spkr_radius^2-y^2);
    th_c = atan2(y,-x_);
    th = th_c;
    spkr_spacing = []; %Auto-calculate spacing
elseif strcmp(array_type, 'line')
    x_=spkr_radius;
    th_c = 180;
    th = atan2(-0.6,-spkr_radius);
    spkr_spacing = 0.001; %1mm spacing between adjacent loudspeakers
end
other = { ...
    'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
    'reproduction_radius',          1.0, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 0, ...
    'unattended_weight',            0.05, ...
    'brightzone_radius',            0.3, ...
    'brightzone_pos_distance',      0.6, ...
    'quietzone_radius',             0.3, ...
    'quietzone_pos_distance',       0.6, ...
    'maximum_frequency',            8000, ...
    'angleof_loudspeakerarrcentre', 180, ...
    'loudspeaker_object',           Parametric_Synthesis.parametric_soundfield };
main_layout = { ...
    'brightzone_pos_angle',        90, ...
    'quietzone_pos_angle',         -90, ...
    'brightzone_source_angle',     0, ...
    'brightzone_source_type',      'pw'};
masker_layout = { ...
    'brightzone_pos_angle',        -90, ...
    'quietzone_pos_angle',         90, ...
    'brightzone_source_angle',     180, ...
    'brightzone_source_type',      'ps'};
loudspeaker_layout = { ...
    'angleto_firstloudspeaker',      90, ...
    'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
    'numberof_loudspeakers',         N_spkrs, ...
    'loudspeaker_model',             'Genelec 8010A', ...
    'loudspeaker_radius',            spkr_radius, ...
    'loudspeaker_spacing',           spkr_spacing, ...
    'speaker_array_type',            array_type, ...
    'brightzone_source_dist',        x_, ...
    other{:}};
main_loudspeaker_layout = loudspeaker_layout;

    %% Room Setup
    Room_Setup = Room_Acoustics.Room;
    Room_Setup.NoReceivers = 32;
    % % ROOM 1
    % % Anechoic
    Room_Setup.Room_Size = [10 10 10]; %Anechoic
    Room_Setup.Wall_Absorb_Coeff = 1.0;

%%
Main_Setup = Speaker_Setup.createSetup({ main_layout{:}, main_loudspeaker_layout{:}});
Masker_Setup = Speaker_Setup.createSetup({ masker_layout{:}, loudspeaker_layout{:}});

masker_signal_info = signal_info;
masker_signal_info.method = 'ZoneWeightMaskerAliasCtrl';

%%
for noise_mask = [-40 -35 -30 -25 -20 -15 -10 -5 0 ]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    masker_signal_info.L_noise_mask = noise_mask; % dB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Playback from setup, setup_info and system_info and record the result
    % Recordings = Hardware_Control.playMSR_andRec( Main_Setup, Room_Setup, signal_info, system_info );
    %Recordings = Hardware_Control.playMSR_andRec( Main_Setup, Room_Setup, signal_info, system_info, Masker_Setup, masker_signal_info );
    %Hardware_Control.playMSR_andRec( Main_Setup, Room_Setup, signal_info, system_info );
    Hardware_Control.playMSR_andRec( Main_Setup, Room_Setup, signal_info, system_info, Masker_Setup, masker_signal_info );
    
    fprintf('\nFinished noise mask level %d \n\n',noise_mask);
    
end

