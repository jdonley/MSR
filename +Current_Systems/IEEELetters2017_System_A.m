function SR_SYSTEM = IEEELetters2017_System_A()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Main_Setup(1) & Room_Setup(1) => Human talker recorded in room
% Main_Setup(2) & Room_Setup(1) => Loudspeaker wall recorded in room
% Main_Setup(1) & Room_Setup(2) => Human talker recorded by microphone wall
% Main_Setup(2) & Room_Setup(2) => NOT USED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
array_type = '2line';
dipoledist = 343/(2*pi*2000);
spkr_radius = 1.5;% - dipoledist/2;
N_spkrs = 24 * 2; % Times 2 for dipole

geometry = 'rectangular';
% geometry = 'circle';

dimensions = 3;

%% Room Geometry
Room_Setup = Room_Acoustics.Room;

% % ROOM 1
% % Anechoic
% Room_Setup = Room_Setup.setRoomSize( [10 10 10] ); %Anechoic
% Room_Setup = Room_Setup.setRoomSize( [4 9 3] ); % 35.G46e
%Room_Setup = Room_Setup.setRoomSize( [8 10 3] ); % 6.107
%Room_Setup = Room_Setup.setRoomSize( [9 14 3] ); % Out to lunch (Cafe)
Room_Setup = Room_Setup.setRoomSize( [3 3 3] ); % Custom [y, x, z]

Room_Setup(1).SystemType = 'Transmit';
Room_Setup(1) = Room_Setup(1).setReproductionCentre( ...
    Room_Setup(1).Room_Size/2 ); % Centre of room [y, x, z]

Room_Setup(1) = Room_Setup(1).setWall_Absorb_Coeff( ...
    [0.2, [1, 1, 1, 1, 1]*1.0]);
% Room_Setup = Room_Setup.setWall_Absorb_Coeff(1.0); % Anechoic (for testing) (comment out otherwise)

% Room_Setup.Reflection_Order = 1; % Order of image sources to compute

Room_Setup(1).NoReceivers = 1; % Number per zone
% If positions are not specified then the positions will be randomised
Room_Setup(1).ReceiverPositions = ...
    Room_Setup(1).Room_Size([2 1 3])/2 + [-0.5 0 0];

%% Multizone Soundfield Geometry and Loudspeaker Array
By = 0.0;
Bx = 0.0;
Qy = 0.0;
Qx = 0.0;




spkrLen    = Room_Setup(1).Room_Size(1);
BZr = spkrLen/2;
QZr = spkrLen/2;
ReproSize  = [spkrLen spkrLen];  % [Width, Height] 
BrightSize = ReproSize;          % [Width, Height] 
QuietSize  = ReproSize;          % [Width, Height] 

srcX = 0.0;
srcY = 0.0;
[srcA,srcD] = cart2pol(srcX,srcY);

imgsrcX = -Room_Setup(1).Room_Size(2) - srcX;
imgsrcY = srcY;
[imgsrcA,imgsrcD] = cart2pol(imgsrcX,imgsrcY);


[Ba,Br] = cart2pol(Bx,By);
[Qa,Qr] = cart2pol(Qx,Qy);

field_layout = { ...
    'brightzone_pos_angle',        Ba/pi*180, ...
    'quietzone_pos_angle',         Qa/pi*180, ...
    'brightzone_source_angle',     imgsrcA/pi*180, ... -atand(0.6/1.3), ...
    'brightzone_source_dist',      imgsrcD, ...
    'brightzone_source_type',      'ps', ...
    'brightzone_geometry',         geometry, ...
    'dimensionality',              dimensions, ...
    'brightzone_size',             BrightSize, ...
    'quietzone_geometry',          geometry, ...
    'quietzone_size',              QuietSize, ...
    'room_size',                   Room_Setup.Room_Size, ...
    'reproduction_geometry',       geometry, ...
    'reproduction_size',           ReproSize};

if strcmpi(array_type, 'circle')
    spkr_spacing = []; %Auto-calculate spacing
elseif contains(lower(array_type), 'line')
    spkr_spacing = 0.001; %1mm spacing between adjacent loudspeakers
end

    loudspeaker_layout = { ...
        'angleto_firstloudspeaker',      90, ...
        'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
        'numberof_loudspeakers',         N_spkrs, ...
        'loudspeaker_model',             'Genelec 8010A', ...
        'loudspeaker_radius',            spkr_radius, ...
        'loudspeaker_spacing',           spkr_spacing, ...
        'speaker_array_type',            array_type, ...
        'angleof_loudspeakerarrcentre', 180, ...
        'quiet_weight',                 1e2};


N_tlkrs = 1;
talker_layout = { ...
        'angleto_firstloudspeaker',      srcA/pi*180, ...
        'angleof_loudspeakerarc',        0 , ...
        'numberof_loudspeakers',         N_tlkrs, ...
        'loudspeaker_model',             'Human', ...
        'loudspeaker_radius',            srcD, ...
        'loudspeaker_spacing',           0, ...
        'speaker_array_type',            'line', ...
        'angleof_loudspeakerarrcentre', srcA/pi*180, ...
        'quiet_weight',                 1e2};
    
% imagesrc_layout = { ...
%         'angleto_firstloudspeaker',      imgsrcA/pi*180, ...
%         'angleof_loudspeakerarc',        0 , ...
%         'numberof_loudspeakers',         N_tlkrs, ...
%         'loudspeaker_model',             'Human', ...
%         'loudspeaker_radius',            imgsrcD, ...
%         'loudspeaker_spacing',           0, ...
%         'speaker_array_type',            'line', ...
%         'angleof_loudspeakerarrcentre', imgsrcA/pi*180, ...
%         'quiet_weight',                 1e2};

c = 343;
Falias = Broadband_Tools.getAliasingFrequency(...
    Speaker_Setup.createSetup({ ...
    field_layout{:}, ...
    loudspeaker_layout{:}, ...
    'brightzone_radius',            BZr, ...
    'brightzone_pos_distance',      Br, ...
    'quietzone_radius',             QZr, ...
    'quietzone_pos_distance',       Qr})) * c / (2*pi);
% Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem.
% We choose 100 for good measure unless memory requirements are an issue due to large room reproductions.
% Taking into consideration the aliasing frequency could significantly improve LUT building performance as well.
field_res = 50;   
field_res = ceil(2*Falias/c / 10)*10; % Round up to nearest multiple of ten


% The talker (speaker) setup
Main_Setup(1) = Speaker_Setup.createSetup({...
    'frequency',                    1000, ...
    field_layout{:}, ...
    talker_layout{:}, ...
    'resolution',                   field_res, ...
    'reproduction_radius',          BZr, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 0, ...
    'unattended_weight',            0.05, ...
    'brightzone_radius',            BZr, ...
    'brightzone_pos_distance',      Br, ...
    'quietzone_radius',             QZr, ...
    'quietzone_pos_distance',       Qr, ...
    'maximum_frequency',            8000});

% The cancellation loudspeaker setup
Main_Setup(2) = Speaker_Setup.createSetup({...
    'frequency',                    1000, ...
    field_layout{:}, ...
    loudspeaker_layout{:}, ...
    'resolution',                   field_res, ...
    'reproduction_radius',          BZr, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 0, ...
    'unattended_weight',            0.05, ...
    'brightzone_radius',            BZr, ...
    'brightzone_pos_distance',      Br, ...
    'quietzone_radius',             QZr, ...
    'quietzone_pos_distance',       Qr, ...
    'maximum_frequency',            8000 });


Falias = 2000;
Main_Setup(2).DipoleDistance =  c / (2*pi*Falias);


Main_Setup(2).ExtendedField = false;
Main_Setup(2).Origin   = Room_Setup.Reproduction_Centre(1:2) - Room_Setup.Room_Size(1:2)./2;

Masker_Setup = [];

%% Receiving Microphone Room Setup
Room_Setup(2) = Room_Setup;
Room_Setup(2).SystemType = 'Receive';

M = Main_Setup(2).Loudspeaker_Count/2; % Number per zone
Room_Setup(2).NoReceivers = M;

spkrLocs = Main_Setup(2).Loudspeaker_Locations(1:M,:);
[spkrLocs(:,1),spkrLocs(:,2)]=pol2cart(spkrLocs(:,1),spkrLocs(:,2));

Room_Setup(2).ReceiverPositions = ...
    [spkrLocs, ...
    zeros(Room_Setup(2).NoReceivers,1)] ...
    + Room_Setup(2).Reproduction_Centre([2 1 3]) ...
    + [dipoledist/2 0 0];

Room_Setup(2) = Room_Setup(2).setReceiverDirectivity('cardioid');
Room_Setup(2).ReceiverOrientations = zeros(M,2);



%% Signal Setup and Path Info
signal_info.c = c; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = ...
                    ...4 ...
                    ...8 ...
                    ...12 ...
                    ...16 ...
                    ...20 ...
                    ...24 ...
                    ...28 ...
                    32 ...
                    * 1e-3 * signal_info.Fs;%1024;% Number of fft components
signal_info.Nfft_optimal = signal_info.Nfft; % The optimal Nfft length for soundfield suppression
signal_info.time_delay = [0] *1e-3; % Seconds %If empty the time delay will based on the frame length
signal_info.system_time_delay = [0] *1e-3; % Seconds %The artificial system time delay between receiving a signal and reproducing the soundfield
signal_info.magnitudeADJ = 0.8;
signal_info.overlap = 0.50;
signal_info.zeropadtime = signal_info.Nfft / signal_info.Fs; % miliseconds of zero padding (i.e. maximum time shift before circular convolution)
signal_info.predict_buff = 2; % Length of prediction buffer as a percentage of a full frame prediction (i.e. 10 is %1000 of frame length)
signal_info.AR_method = 'lpc';
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.f_low_meas = 100; % Hz %Minimum loudspeaker response
signal_info.f_high_meas = 7000; % Hz %Maximum frequency with accurate response at given sampling rate
signal_info.L_noise_mask = [-inf]; % dB
signal_info.recording_type = {'simulated'}; % The type of recordings to be analysed
signal_info.weight = 1; % This can be auto-calculated for maximum contrast by setting to 'Auto'
signal_info.method = ''; % Default empty (temporary variable)
signal_info.methods_list ... % List of methods to synthesize
    = {'Clean'; ... % There should be at least two entries here
       'BoundaryCancel'};
   
signal_info.methods_list_clean = [1;2]; %Indices from the "methods_list" of the clean signals
signal_info.methods_list_masker = [0]; %Indices from the "methods_list" of the maskers, different hybrids are separated by columns
% ( e.g. [2,3;4,0;6,7] is two hybrids, the first is 2&4&6, the second is 3&7, indices < 1 are ignored)
signal_info.reference = false; % False or True to record reference signal
signal_info.reference_channel = 1; %Some arbitrary reference signal channel
signal_info.rir_duration = 0.5; % Room Impulse Response length in seconds
signal_info.input_filename = [];
% signal_info.speech_filepath = '+Miscellaneous\+Speech_Files\';
signal_info.speech_filepath = '+Miscellaneous\+TestAudio_Files\';

% signal_info.speech_filepath = '+Miscellaneous\+Speech_File_Test\';
%signal_info.speech_filepath = '+Miscellaneous\+Noise_Files\';
%signal_info.speech_filepath = '+Miscellaneous\+STIPA_Test\';
%signal_info.speech_filepath = '+Miscellaneous\+Impulse_Response\';
%signal_info.speech_filepath = '+Miscellaneous\+Sine_Sweep\';

signal_info.InverseFilter_filepath = '+Miscellaneous\+TestAudio_Files_InvFilts\';

%% System Setup
system_info.dev_model = 'ASIO Hammerfall DSP';
system_info.fs = 48000;
system_info.f_low = 100; % Hz %Minimum calibration frequency
system_info.f_high = 10000; % Hz %Maximum calibration frequency
system_info.playbackChannels = ...
    [ 1  2  3  4  5  6  7  8 ...
    9 10 11 12 13 14 15 16 ...
    17 18];% 19 20 21 22 23 24];

% First half recordings in Bright zone, second half in Quiet zone
system_info.recordChannels = ...
    [ 1 2 ];

system_info.calibrationRecChannel = ...
    [ 3 ];

system_info.Sweep_Length = 10; %Seconds
system_info.Sweep_EndBuffers = 1; %Seconds
system_info.Calibration_FiltLen = 0.5; %Seconds %Filter length
system_info.Calibration_FiltReg = [60 -6]; %[passband_gain, stopband_gain] (dB)

system_info.sc = '_'; % Separating character for ascii paths
system_info.Drive = ['Z:' filesep]; % Database drive (storage drive)
system_info.Filter_dir = ['+Speaker_Signals' filesep 'CalibratedSweeps' filesep];
system_info.FilterData_dir = ['+Calibration' system_info.sc 'Data' filesep '+Filters' filesep];
system_info.CalibrationRec_dir = ['+Calibration' system_info.sc 'Data' filesep '+Recordings' filesep];
system_info.Calibrated_Signals_dir = ['+Calibrated' system_info.sc 'Speaker_Signals' filesep];

system_info.LUT_frequencies = (signal_info.Nfft + signal_info.zeropadtime * signal_info.Fs)/2;
% system_info.LUT_weights = 32;
system_info.LUT_weights = 1; % Number of Look-Up table multizone soundfield weights
% system_info.LUT_weight_range = [1e-2 1e4]; % [Minimum Maximum] LUT weight
system_info.LUT_weight_range = 1; % [Minimum Maximum] or a single LUT weight

system_info.LUT_resolution = [num2str(system_info.LUT_frequencies) 'f' ...
                              system_info.sc ...
                              num2str(system_info.LUT_weights) 'w'];
                          
%% Analysis Information
analysis_info.Measures    = {'RIR'};
analysis_info.Result_Type = {'Predicted Signal'; ...
                             'Actual Signal'};
analysis_info.SignalNames = {'Active Wall Off'; ...
                             'Active Wall On'};
analysis_info.Nfft = 32 * 1e-3 * signal_info.Fs;%1024;% Number of fft components
analysis_info.f_low = 150; % Hz
analysis_info.f_high = 2000; % Hz
% analysis_info.f_high = Falias; % Hz


%% Publication Figure Setup Information
publication_info.DocumentPath = 'tex\latex\IEEE_Letters2017';


publication_info.FigureName = 'IEEE_Letters2017_A';
publication_info.axes_NumTicks = [8 5]; % Number of ticks [NumXticks NumYticks];
publication_info.axes_limitBufs = [0.00 0.05]; % axis limits buffer in percentage [width, height]
publication_info.axes_Scales = {'lin'; 'lin'}; % Axis scales
% publication_info.XTicks_override = [0.156, 1, 8];
publication_info.leg_MarkerSize = 8;

publication_info.FigureTitle = '';%'Active Speech Attenuation';
publication_info.print_fmt = 'pdf'; %figure image file format
publication_info.print_res = 600; %rastered graphic DPI

publication_info.figure_width = 88.9/10;% + 6.35/10 + 88.9/10; %Figure width in centimeters %IEEE full text width
publication_info.figure_aspect_ratio = 6/3; %Full figure aspect ration width/height
publication_info.axis_aspect_ratio = [1 0.4]; %Single axis aspect ration [width height]
publication_info.axes_gap = [0.5 0.5]; %Gap between axes [gap_height gap_width] %centimeters
publication_info.axes_margins_height = [1 1]; %Axes height margins [lower upper]  %centimeters
publication_info.axes_margins_width = [1 1]; %Axes width margins [left right]  %centimeters
publication_info.axes_grid = 'off'; % Show axes grid ('on', 'minor' or 'off')
publication_info.axes_gridMinorY = 'off'; % Override minor grid for Y axis ('on' or 'off')
publication_info.axes_gridMinorX = 'off'; % Override minor grid for X axis ('on' or 'off')
publication_info.sigRounding = 3; % number of significant figures rounding
publication_info.axes_tickdir = 'both'; % Axes tick direction(s) ('in', 'out' or 'both')

publication_info.FontSize = 9;  % Font size of text in figure
publication_info.FontName = 'Times'; % Font name of text in figure
publication_info.NumbersFontName = 'fixedwidth'; % Font name of numbers in figure
publication_info.Interpreter = 'latex'; % Interpreter of text in figure
publication_info.LaTeX_FontFamily = 'cmr'; % Font name of text in figure
publication_info.LaTeX_NumbersFontFamily = 'cmtt'; % Font name of numbers in figure
publication_info.lineWid = 0.5; % PDF line widths

%% Clear all except variables needed
% Create soundfield reproduction system structure
SR_SYSTEM = Current_Systems.SR_System;
SR_SYSTEM.Room_Setup  = Room_Setup;
SR_SYSTEM.Main_Setup  = Main_Setup;
SR_SYSTEM.Masker_Setup  = Masker_Setup;
SR_SYSTEM.signal_info = signal_info;
SR_SYSTEM.system_info = system_info;
SR_SYSTEM.analysis_info = analysis_info;
SR_SYSTEM.publication_info = publication_info;

end
