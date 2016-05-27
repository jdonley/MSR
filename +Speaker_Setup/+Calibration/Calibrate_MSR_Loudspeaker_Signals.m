clc;clear;


%% Setup Information
%Flip loudspeaker order (effectively flips entire setup) (false if not needed)
MirrorSetup = false;

system_info.sc = '_'; % Separating character for ascii paths
system_info.Drive = 'Z:\';
Calibration_Path = [system_info.Drive '+Calibration_Data' filesep];
Filters_Dir = ['+Filters' filesep];
Filter_path = [Calibration_Path Filters_Dir];
Calibrated_Signals_dir = '+Calibrated_Speaker_Signals\';


LUT_resolution = '512f_32w';

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


array_type = 'circle';
spkr_radius = 1.3;
N_spkrs = 24;

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
loudspeaker_layout = [{ ...
    'angleto_firstloudspeaker',      90, ...
    'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
    'numberof_loudspeakers',         N_spkrs, ...
    'loudspeaker_model',             'Genelec 8010A', ...
    'loudspeaker_radius',            spkr_radius, ...
    'loudspeaker_spacing',           spkr_spacing, ...
    'speaker_array_type',            array_type, ...
    'brightzone_source_dist',        x_, ...
    },other];
main_loudspeaker_layout = loudspeaker_layout;

Main_Setup = Speaker_Setup.createSetup([ main_layout, main_loudspeaker_layout]);



Setup = Speaker_Setup.createSetup([ main_layout, loudspeaker_layout]);
%%
noise_levels = signal_info.L_noise_mask;

% Comment out the following three lines to calibrate clean speech files
Setup = Speaker_Setup.createSetup({ masker_layout{:}, loudspeaker_layout{:}});
signal_info.method = 'ZoneWeightMaskerAliasCtrl';
noise_levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20];

for noise_mask = noise_levels
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal_info.L_noise_mask = noise_mask; % dB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Get Paths and files
    [Spkr_path,~,~,~,~,~,path_ext] = ...
        Broadband_Tools.getLoudspeakerSignalPath( Setup, signal_info, LUT_resolution, system_info.Drive, 'new');
    
    files = Tools.getAllFiles(Spkr_path);
    files = sort(files);
    
    filter_location = Tools.getAllFiles(Filter_path);
    filter_location = sort(filter_location);
    
    
    %% Load Fiters
    filts = load( filter_location{1} ); % 1st element should be the newest (most recent) set of filters
    
    if MirrorSetup
        filts.EQ = flip(filts.EQ,2);
    end
    
    %%
    fileName_prev = '';
    Speaker_Signals = [];
    
    F=length(files);
    for f = 1:F
        
        [~, fileName, fileExt] = fileparts(files{f});
        if isempty( strfind( fileName, 'Original' ) ) % If not an Original audio file
            try
                y = audioread(files{f});
            catch err
                if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                    continue; % Skip unsupported files
                end
            end
            
            % Get the file number and file name
            [fnumflip,fnameflip] = strtok( flip(fileName), system_info.sc );
            Current_Spkr_Num = str2double( flip( fnumflip ) );
            fileName_curr = flip(sscanf(fnameflip,['%*[' system_info.sc ']%s']));
            
            if ~(isempty(fileName_prev) || strcmp( fileName_curr, fileName_prev))
                Speaker_Signals = [];
            end
            
            Speaker_Signals(Current_Spkr_Num,:) = y;
            
            if sum( any( Speaker_Signals, 2 ) ) == Main_Setup.Loudspeaker_Count % If we have a complete group of speaker signals
                Speaker_Signals = Speaker_Signals';
                Speaker_Signals = flip(Speaker_Signals,2); % Speakers are in reverse order
                
                %% Read original file
                if isempty(strfind(signal_info.method, 'Masker')) % If the signals are not Maskers then read the Original audio
                    try
                        orig = audioread( [Spkr_path, fileName_curr, system_info.sc, 'Original', fileExt] ); % Assumes same file extension as loudspeaker audio files
                    catch err
                        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                            error('Error reading the Original audio file for the reproduction.'); % Skip unsupported files
                        end
                    end
                else
                    orig=[];
                end
                
                %% Upsample to the filter sampling frequency
                up_factor = filts.fs/signal_info.Fs;
                [b,a] = butter( 9, 1/up_factor );
                loudspeaker_signals_upsampled = zeros(size(Speaker_Signals,1)*up_factor, Main_Setup.Loudspeaker_Count);
                for spkr = 1:Main_Setup.Loudspeaker_Count
                    loudspeaker_signals_upsampled(:,spkr) = ...
                        filter( b, a, ...
                        interp( Speaker_Signals(:,spkr), up_factor ) );
                end
                
                if ~isempty(orig)
                    orig_upsampled = filter( b, a, interp( orig, up_factor ) );
                end
                
                
                %% Calibrate
                loudspeaker_signals_calibrated = Tools.fconv( ...
                    loudspeaker_signals_upsampled, ...
                    filts.EQ );                
                if ~isempty(orig)
                    orig_calibrated = Tools.fconv( ...
                        orig_upsampled, ...
                        [1; zeros(size(filts.EQ,1)-1,1)] );
                end
                
                %% Save
                % Normalise
                scale_value = 1 ./ max(abs(loudspeaker_signals_calibrated(:)));
                loudspeaker_signals_calibrated = loudspeaker_signals_calibrated .* scale_value;
                if ~isempty(orig)
                    orig_calibrated = orig_calibrated ./ max(abs(orig_calibrated(:)));
                end
                
                output_dir = [system_info.Drive Calibrated_Signals_dir path_ext];
                
                if ~exist(output_dir,'dir'); mkdir(output_dir); end
                
                for spkr = 1:Main_Setup.Loudspeaker_Count
                    audiowrite([output_dir fileName_curr system_info.sc 'Upsampled' system_info.sc num2str(spkr) fileExt], loudspeaker_signals_calibrated(:,spkr), filts.fs);
                end
                
                audiowrite([output_dir fileName_curr system_info.sc 'Upsampled' fileExt], loudspeaker_signals_calibrated, filts.fs);
                save([output_dir fileName_curr system_info.sc 'Upsampled' '.mat'], 'scale_value');
                         
                if ~isempty(orig)
                    audiowrite([output_dir fileName_curr system_info.sc 'Original' fileExt], orig_calibrated, filts.fs );
                end
                
            end
            
            fileName_prev = fileName_curr;
        end
    end
    
    
end







