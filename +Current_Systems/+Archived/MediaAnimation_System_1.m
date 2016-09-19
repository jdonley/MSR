function SR_SYSTEM = MediaAnimation_System_1()

array_type_ = {'line';'circle'};
N_spkrs_ = 22;

for at = 1:2
    array_type = array_type_{at};
    for Nsp = 2:N_spkrs_
        N_spkrs = Nsp;
        
        spkr_type  = 'Dynamic';
        spkr_radius = 1.3;
        
        %% Room Geometry
        Room_Setup = Room_Acoustics.Room;
        Room_Setup.NoReceivers = 32;
        % % ROOM 1
        % % Anechoic
        Room_Setup.setRoomSize( [10 10 10] ); %Anechoic
        %Room_Setup.setRoomSize( [4 9 3] ); % 35.G46e
        %Room_Setup.setRoomSize( [8 10 3] ); % 6.107
        %Room_Setup.setRoomSize( [9 14 3] ); % Out to lunch (Cafe)
        
        Room_Setup.setReproductionCentre( Room_Setup.Room_Size .* [0.5 0.5 0.5] ); % Centre of room
        
        Room_Setup = Room_Setup.setWall_Absorb_Coeff(1.0);
        
        %% Multizone Soundfield Geometry and Loudspeaker Array
        gemoetrical_layout = { ...
            'brightzone_pos_angle',        90, ...
            'quietzone_pos_angle',         -90, ...
            'brightzone_source_angle',     0, ... -atand(0.6/1.3), ...
            'brightzone_source_type',      'pw'};
        masker_layout = { ...
            'brightzone_pos_angle',        -90, ...
            'quietzone_pos_angle',         90, ...
            'brightzone_source_angle',     0, ... -atand(0.6/1.3), ...
            'brightzone_source_type',      'pw'};
        Para_Spkr = Parametric_Synthesis.parametric_soundfield;
        Para_Spkr.P1 = db2mag( 100 ); % 100dB amplitude parametric array loudspeaker
        Para_Spkr.P2 = db2mag( 100 ); % 100dB secondary amplitude
        if strcmpi(array_type, 'Circle')
            [x,y] = pol2cart(-90/180*pi, 0.6);
            x_ = sqrt(spkr_radius^2-y^2);
            th_c = atan2(y,-x_);
            th = th_c;
            spkr_spacing = []; %Auto-calculate spacing
            angleof_loudspeakerarc = 180 * N_spkrs/(N_spkrs-1);
        elseif strcmpi(array_type, 'Line')
            x_=spkr_radius;
            th_c = 180;
            th = atan2(-0.6,-spkr_radius);
            spkr_spacing = (spkr_radius*2-N_spkrs*0.121)/(N_spkrs-1); %1mm spacing between adjacent loudspeakers
            angleof_loudspeakerarc = 218.6;
        end
        if strcmpi(spkr_type, 'Dynamic')
            loudspeaker_layout = { ...
                'angleto_firstloudspeaker',      90, ...
                'angleof_loudspeakerarc',        angleof_loudspeakerarc , ...
                'numberof_loudspeakers',         N_spkrs, ...
                'loudspeaker_model',             'Genelec 8010A', ...
                'loudspeaker_radius',            spkr_radius, ...
                'loudspeaker_spacing',           spkr_spacing, ...
                'speaker_array_type',            array_type, ...
                'brightzone_source_dist',        x_, ...
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
        Main_Setup(Nsp-1,at) = Speaker_Setup.createSetup({...
            'frequency',                    1000, ...
            gemoetrical_layout{:}, ...
            loudspeaker_layout{:}, ...
            'resolution',                   200, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
            'reproduction_radius',          1.0, ...
            'bright_weight',                1.0, ...
            'unattended_weight',            0.05, ...
            'brightzone_radius',            0.3, ...
            'brightzone_pos_distance',      0.6, ...
            'quietzone_radius',             0.3, ...
            'quietzone_pos_distance',       0.6, ...
            'maximum_frequency',            8000, ...
            'loudspeaker_object',           Para_Spkr });
        Masker_Setup = Speaker_Setup.createSetup({...
            'frequency',                    1000, ...
            masker_layout{:}, ...
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
            'loudspeaker_object',           Para_Spkr });
        
    end
end

%% Signal Setup and Path Info
signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.f_low_meas = 100; % Hz %Minimum loudspeaker response
signal_info.f_high_meas = 7000; % Hz %Maximum frequency with accurate response at given sampling rate
signal_info.L_noise_mask = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20]; % dB
signal_info.recording_type = {'simulated'; 'realworld'}; % The type of recordings to be analysed
signal_info.weight = 100; % This can be auto-calculated for maximum contrast by setting to 'Auto'
signal_info.method = ''; % Default empty (temporary variable)
signal_info.methods_list ... % List of methods to synthesize
    = {'NoMask'};%, ... % Speech Signal
%'FlatMasker', ...
%'ZoneWeightMaskerAliasCtrl'};
signal_info.methods_list_clean = [1]; %Indices of the clean signals
signal_info.methods_list_masker = [];%[2, 3]; %Indices of the maskers, different hybrids are separated by columns
% ( e.g. [2,3;4,0;6,7] is two hybrids, the first is 2&4&6, the second is 3&7, indices < 1 are ignored)
signal_info.reference = false; % False or True to record reference signal
signal_info.reference_channel = 1; %Some arbitrary reference signal channel
signal_info.rir_duration = 0.5; % Room Impulse Response length in seconds
signal_info.input_filename = [];
signal_info.speech_filepath = '+Miscellaneous\+Speech_Files\';
%signal_info.speech_filepath = '+Miscellaneous\+Speech_File_Test\';
%signal_info.speech_filepath = '+Miscellaneous\+Noise_Files\';
%signal_info.speech_filepath = '+Miscellaneous\+TestAudio_Files\';
%signal_info.speech_filepath = '+Miscellaneous\+STIPA_Test\';
%signal_info.speech_filepath = '+Miscellaneous\+Impulse_Response\';
%signal_info.speech_filepath = '+Miscellaneous\+Sine_Sweep\';

%% System Setup
system_info.dev_model = 'ASIO Hammerfall DSP';
system_info.fs = 48000;
system_info.f_low = 100; % Hz %Minimum calibration frequency
system_info.f_high = 10000; % Hz %Maximum calibration frequency
system_info.playbackChannels = ...
    [ 1  2  3  4  5  6  7  8 ...
    9 10 11 12 13 14 15 16 ...
    17 18 19 20 21 22 23 24];

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

system_info.LUT_frequencies = 512;
% system_info.LUT_weights = 32;
system_info.LUT_weights = 2; % First weight is always zero
% system_info.LUT_weight_range = [1e-2 1e4]; % [Minimum Maximum] LUT weight
system_info.LUT_weight_range = 1e2; % [Minimum Maximum] or a single LUT weight

system_info.LUT_resolution = [num2str(system_info.LUT_frequencies) 'f' ...
    system_info.sc ...
    num2str(system_info.LUT_weights) 'w'];

%% Publication Figure Setup Information
publication_info.DocumentPath = 'tex\latex\IEEE_Trans2016';
publication_info.FigureName = 'IEEE_Trans2016_A';
publication_info.FigureTitle = 'Quality and Intelligibility Results';
publication_info.print_fmt = 'pdf'; %figure image file format
publication_info.print_res = 600; %rastered graphic DPI

publication_info.figure_width = 88.9/10 + 6.35/10 + 88.9/10; %Figure width in centimeters %IEEE full text width
publication_info.figure_aspect_ratio = 6/3; %Full figure aspect ration width/height
publication_info.axis_aspect_ratio = [1 0.4]; %Single axis aspect ration [width height]
publication_info.axes_gap = [0.5 0.5]; %Gap between axes [gap_height gap_width] %centimeters
publication_info.axes_margins_height = [1 1]; %Axes height margins [lower upper]  %centimeters
publication_info.axes_margins_width = [1 1]; %Axes width margins [left right]  %centimeters
publication_info.axes_grid = 'minor'; % Show axes grid ('on', 'minor' or 'off')
publication_info.axes_NumTicks = [5 6]; % Number of ticks [NumXticks NumYticks];
publication_info.axes_limitBufs = [0.02 0.02]; % axis limits buffer in percentage [width, height]
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
SR_SYSTEM = Current_Systems.SR_System;
SR_SYSTEM.Room_Setup  = Room_Setup;
SR_SYSTEM.Main_Setup  = Main_Setup;
SR_SYSTEM.Masker_Setup  = Masker_Setup;
SR_SYSTEM.signal_info = signal_info;
SR_SYSTEM.system_info = system_info;
SR_SYSTEM.publication_info = publication_info;

end
