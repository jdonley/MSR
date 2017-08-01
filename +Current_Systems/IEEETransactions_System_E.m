function SR_SYSTEM = IEEETransactions_System_E()


array_type = {'circle','line'}; % 'circle' or 'line'

SourceAngleSetIndice = 2; % 1, 2 or 3

N_spkrs = 24; % 16, 24, 32 or 149

% lambda_g = [0.0, 0.5, 1.0]; % between 0 and 1
% lambda_gTxt = num2str(lambda_g.','%0.2f');
% ZWMAC = 'ZoneWeightMaskerAliasCtrl';
% masker_types = [{'FlatMasker'}, ...
%     mat2cell( ...
%     reshape([repmat(ZWMAC,numel(lambda_g),1),repmat(lambda_gTxt,1,1)].',1,[]), ...
%     1,(numel(ZWMAC)+size(lambda_gTxt,2))*ones(numel(lambda_g),1))];

rec_types = {'simulated'; 'realworld'}; % The type of recordings to be analysed

spkr_type  = 'Dynamic';
spkr_radius = 1.3;
% spkr_radius = 1.32;
% spkr_radius = 1.34;
% spkr_radius = 1.36;

dimensions = 3;

%% Room Geometry
Room_Setup = Room_Acoustics.Room;
Room_Setup.NoReceivers = 32;
% % ROOM 1
% % Anechoic
Room_Setup = Room_Setup.setRoomSize( [10 10 10] ); %Anechoic
%Room_Setup = Room_Setup.setRoomSize( [4 9 3] ); % 35.G46e
%Room_Setup = Room_Setup.setRoomSize( [8 10 3] ); % 6.107
%Room_Setup = Room_Setup.setRoomSize( [9 14 3] ); % Out to lunch (Cafe)

Room_Setup = Room_Setup.setReproductionCentre( Room_Setup.Room_Size .* [0.5 0.5 0.5] ); % Centre of room

Room_Setup = Room_Setup.setWall_Absorb_Coeff(1.0);

%% Multizone Soundfield Geometry and Loudspeaker Array
% N_sets = numel(lambda_g)*numel(array_type);
Bx =  0;
By =  0.6;
Qx =  0;
Qy = -0.6;

for arr = 1:numel(array_type)
    array_type_ = array_type{arr};
    
    I = arr;
    
    if strcmpi(array_type_,'circle')
        switch SourceAngleSetIndice
            case 1
                Theta    =  0;
                Vartheta = -90+acosd( (abs(By)+abs(Qy)) / sqrt(abs(2*By*Qy)+Qy^2+spkr_radius^2) );
            case 2
                Theta    =  atand(mean(abs([By,Qy]))/spkr_radius);
                Vartheta = -atand(mean(abs([By,Qy]))/spkr_radius);
            case 3
                Theta    =  90-acosd( (abs(By)+abs(Qy)) / sqrt(abs(2*By*Qy)+By^2+spkr_radius^2) );
                Vartheta =  0;
        end
    elseif strcmpi(array_type_,'line')
        switch SourceAngleSetIndice
            case 1
                Theta    =  0;
                Vartheta = -90+atand( spkr_radius / (abs(By)+abs(Qy)) );
            case 2
                Theta    =  atand(mean(abs([By,Qy]))/spkr_radius);
                Vartheta = -atand(mean(abs([By,Qy]))/spkr_radius);
            case 3
                Theta    =  90-atand( spkr_radius / (abs(By)+abs(Qy)) );
                Vartheta =  0;
        end
    end
    
    gemoetrical_layout = { ...
        'brightzone_pos_angle',        90, ...
        'quietzone_pos_angle',         -90, ...
        'brightzone_source_angle',     Theta, ...
        'brightzone_source_dist',      sqrt(0.6^2+1.3^2), ...
        'brightzone_source_type',      'pw'};
    masker_layout = { ...
        'brightzone_pos_angle',        -90, ...
        'quietzone_pos_angle',         90, ...
        'brightzone_source_angle',     Vartheta, ...
        'brightzone_source_dist',      sqrt(0.6^2+1.3^2), ...
        'brightzone_source_type',      'pw'};
    Para_Spkr = Parametric_Synthesis.parametric_soundfield;
    Para_Spkr.P1 = db2mag( 100 ); % 100dB amplitude parametric array loudspeaker
    Para_Spkr.P2 = db2mag( 100 ); % 100dB secondary amplitude
    if strcmpi(array_type_, 'Circle')
        [x,y] = pol2cart(-90/180*pi, 0.6);
        x_ = sqrt(spkr_radius^2-y^2);
        th_c = atan2(y,-x_);
        th = th_c;
        spkr_spacing = []; %Auto-calculate spacing
    elseif strcmpi(array_type_, 'Line')
        x_=spkr_radius;
        th_c = 180;
        th = atan2(-0.6,-spkr_radius);
        spkr_spacing = 0.001; %1mm spacing between adjacent loudspeakers
    end
    if strcmpi(spkr_type, 'Dynamic')
        loudspeaker_layout = { ...
            'angleto_firstloudspeaker',      90, ...
            'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
            'numberof_loudspeakers',         N_spkrs, ...
            'loudspeaker_model',             'Genelec 8010A', ...
            'loudspeaker_radius',            spkr_radius, ...
            'loudspeaker_spacing',           spkr_spacing, ...
            'speaker_array_type',            array_type_, ...
            'angleof_loudspeakerarrcentre', 180, ...
            'quiet_weight',                 1e2};
    elseif strcmpi(spkr_type, 'Parametric')
        loudspeaker_layout = {  ...
            'angleto_firstloudspeaker',     atan2d(-0.6,-x_), ...
            'angleof_loudspeakerarrcentre', 180, ... +atand(0.6/1.3), ...
            'loudspeaker_radius',           x_, ... spkr_radius, ...
            'numberof_loudspeakers',        1, ...
            'loudspeaker_model',            'Parametric', ...
            'loudspeaker_spacing',          0.01, ...
            'speaker_array_type',           'line', ...
            'brightzone_source_dist',        x_};
        Para_Spkr = Para_Spkr.set_f1( 40000 );
    end
    Main_Setup(I) = Speaker_Setup.createSetup({...
        'frequency',                    1000, ...
        gemoetrical_layout{:}, ...
        loudspeaker_layout{:}, ...
        'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
        'reproduction_radius',          1.0, ...
        'bright_weight',                1.0, ...
        'unattended_weight',            0.05, ...
        'brightzone_radius',            0.3, ...
        'brightzone_pos_distance',      0.6, ...
        'quietzone_radius',             0.3, ...
        'quietzone_pos_distance',       0.6, ...
        'maximum_frequency',            8000, ...
        'dimensionality',               dimensions, ...
        'loudspeaker_object',           Para_Spkr });
    Masker_Setup = {};
    
end


%% Signal Setup and Path Info
signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.time_delay = 0 *1e-3; % Seconds %If empty the time delay will based on the frame length
signal_info.overlap = 0.5;
signal_info.zeropadtime = 0; % miliseconds of zero padding (i.e. maximum time shift before circular convolution)
signal_info.predict_buff = 0; % Length of prediction buffer as a percentage of a full frame prediction (i.e. 10 is %1000 of frame length)
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.f_low_meas = 100; % Hz %Minimum loudspeaker response
signal_info.f_high_meas = 7000; % Hz %Maximum frequency with accurate response at given sampling rate
signal_info.OctaveBandSpace = 1/12; % Octave band spacing for filtering operations
signal_info.clipLevelAdjust = -35; % dB %RMS level to adjust to in order to avoid clipping (Hope that -35dB RMS level is low enough to avoid clipping upon saving)
signal_info.L_noise_mask = -Inf; % dB
signal_info.recording_type = rec_types; % The type of recordings to be analysed
signal_info.weight = 100; % This can be auto-calculated for maximum contrast by setting to 'Auto'
signal_info.method = ''; % Default empty (temporary variable)


signal_info.methods_list ... % List of methods to synthesize
    = {'NoMask';'NoMask'}; % Speech Signal

signal_info.methods_list_clean = [1 2]; %Indices of the clean signals
signal_info.methods_list_masker = [0];%[2, 3, 4, 5]; %Indices of the maskers, different hybrids are separated by columns
% ( e.g. [2,3;4,0;6,7] is two hybrids, the first is 2&4&6, the second is 3&7, indices < 1 are ignored)
signal_info.reference = false; % False or True to record reference signal


signal_info.reference_channel = 1; %Some arbitrary reference signal channel
signal_info.rir_duration = 0.5; % Room Impulse Response length in seconds
signal_info.UseMeasuredATFs = true; % Use the measured ATFs between loudspeaker and microphone pairs instead of simulated RIRs
signal_info.input_filename = [];
signal_info.inputSignalNorm = true; % Normalise the input signal to RMS value
% signal_info.speech_filepath = '+Miscellaneous\+Speech_Files\';
%signal_info.speech_filepath = '+Miscellaneous\+Speech_File_Test\';
%signal_info.speech_filepath = '+Miscellaneous\+Noise_Files\';
signal_info.speech_filepath = '+Miscellaneous\+TestAudio_Files\';
%signal_info.speech_filepath = '+Miscellaneous\+STIPA_Test\';
%signal_info.speech_filepath = '+Miscellaneous\+Impulse_Response\';
%signal_info.speech_filepath = '+Miscellaneous\+Sine_Sweep\';

signal_info.InverseFilter_filepath = '+Miscellaneous\+TestAudio_Files_InvFilts\';

%% System Setup
system_info.dev_model = 'ASIO Hammerfall DSP';
system_info.fs = 48000;
system_info.f_low = 100; % Hz %Minimum calibration frequency
system_info.f_high = 10000; % Hz %Maximum calibration frequency
system_info.MirrorSetup = false; % Mirrors the loudspeaker setup (flips playback and calibration order)
system_info.playbackChannels = ...
    [ 1  2  3  4  5  6  7  8 ...
    9 10 11 12 13 14 15 16 ...
    17 18 19 20 21 22 23 24];
% system_info.playbackChannels = ...
%     [ 1  ];
% system_info.recordChannels = ...
%     [ 1 2 ];
% system_info.recordChZoneAlloc = ... % Allocation of the recorded channels to their respective zones (1 is for Bright Zone, 2 is for Quiet Zone)
%     [ 1 2 ];

system_info.CurrentSpeakerArrayType = 'line';
% zone = 'bright';
zone = 'quiet';

system_info.recordChannels = ...
    [ 1 2 3 4];
if strcmpi(zone,'bright')
    system_info.recordChZoneAlloc = ... % Allocation of the recorded channels to their respective zones (1 is for Bright Zone, 2 is for Quiet Zone)
        [ 1 1 1 1];
elseif strcmpi(zone,'quiet')
    system_info.MirrorSetup = true;
    system_info.recordChZoneAlloc = ... % Allocation of the recorded channels to their respective zones (1 is for Bright Zone, 2 is for Quiet Zone)
        [ 2 2 2 2];
end



system_info.calibrationRecChannel = ...
    [ 5 ];

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

system_info.LUT_frequencies = 512;
% system_info.LUT_weights = 32;
system_info.LUT_weights = 2; % First weight is always zero
% system_info.LUT_weight_range = [1e-2 1e4]; % [Minimum Maximum] LUT weight
system_info.LUT_weight_range = 1e2; % [Minimum Maximum] or a single LUT weight

system_info.LUT_resolution = [num2str(system_info.LUT_frequencies) 'f' ...
    system_info.sc ...
    num2str(system_info.LUT_weights) 'w'];

%% Analysis Information
analysis_info.Measures = {'SPL'};

analysis_info.SignalNames = {'BZ'; ...
    'QZ'};

analysis_info.Nfft = 32 *1e-3*signal_info.Fs;%1024;% Number of fft components
analysis_info.f_low = 90; % Hz
analysis_info.f_high = 6500; % Hz

%% Publication Figure Setup Information
publication_info.DocumentPath = 'tex\latex\IEEE_Trans2016';
publication_info.FigureName = 'IEEE_Trans2016_E';
publication_info.FigureTitle = 'Zone Sound Pressure Levels - Simulated and Real-World';
publication_info.print_fmt = 'pdf'; %figure image file format
publication_info.print_res = 600; %rastered graphic DPI

publication_info.LatexMacrosFile = 'IEEE_Trans2016_LaTeX_Macros.tex';

publication_info.axes_NumTicks = [8 5]; % Number of ticks [NumXticks NumYticks];
publication_info.axes_limitBufs = [0.00 0.100]; % axis limits buffer in percentage [width, height]
publication_info.axes_Scales = {'log'; 'lin'}; % Axis scales
publication_info.XTicks_override = [0.10, 1, 8];
publication_info.YLim_override = [-40 0]; % Sets the Y limits (range) of the axes
publication_info.leg_MarkerSize = 8;

publication_info.figure_width = 88.9/10;% + 6.35/10 + 88.9/10; %Figure width in centimeters %IEEE full text width
publication_info.figure_aspect_ratio = 8/3; %Full figure aspect ratio width/height
publication_info.axis_aspect_ratio = [1 0.333]; %Single axis aspect ration [width height]
publication_info.axes_gap = [0.5 0.5]; %Gap between axes [gap_height gap_width] %centimeters
publication_info.axes_margins_height = [1 1]; %Axes height margins [lower upper]  %centimeters
publication_info.axes_margins_width = [1 1]; %Axes width margins [left right]  %centimeters
publication_info.axes_grid = 'off'; % Show axes grid ('on', 'minor' or 'off')
publication_info.axes_gridMinorY = 'off'; % Override minor grid for Y axis ('on' or 'off')
publication_info.axes_gridMinorX = 'off'; % Override minor grid for X axis ('on' or 'off')
publication_info.axes_tickdir = 'both'; % Axes tick direction(s) ('in', 'out' or 'both')
publication_info.sigRounding = 3; % number of significant figures rounding

publication_info.FontSize = 9;  % Font size of text in figure
publication_info.FontName = 'Times'; % Font name of text in figure
publication_info.NumbersFontName = 'fixedwidth'; % Font name of numbers in figure
publication_info.Interpreter = 'latex'; % Interpreter of text in figure
publication_info.LaTeX_FontFamily = 'cmr'; % Font name of text in figure
publication_info.LaTeX_NumbersFontFamily = 'cmtt'; % Font name of numbers in figure
publication_info.lineWid = 0.5; % PDF line widths

%% Clear all except variables needed
% Create soundfield reproduction system structure
SR_SYSTEM                   = Current_Systems.SR_System;
SR_SYSTEM.Room_Setup        = Room_Setup;
SR_SYSTEM.Main_Setup        = Main_Setup;
SR_SYSTEM.Masker_Setup      = Masker_Setup;
SR_SYSTEM.signal_info       = signal_info;
SR_SYSTEM.system_info       = system_info;
SR_SYSTEM.analysis_info     = analysis_info;
SR_SYSTEM.publication_info  = publication_info;

end
