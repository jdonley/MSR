function Reverberant_MSR_batchfunc( SYS_or_room_setup, mask_type, main_setup, weight, mask_level, Conv_Type )

SYS_type = 'Current_Systems.SR_System';

%% Initialise
tic;
% Start Parallel Pool
para_pool = parpool;
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%%
if nargin == 1
    if isa(SYS_or_room_setup,SYS_type)
        SYS = SYS_or_room_setup;
        main_setup = SYS.Main_Setup;
        masker_setup = SYS.Masker_Setup;
        room_setup = SYS.Room_Setup;
        signal_info = SYS.signal_info;
        system_info = SYS.system_info;
    else
        error(['Single input argument must be of type: ' SYS_type]);
    end
elseif nargin > 1
    room_setup = SYS_or_room_setup;
    system_info.Drive = 'Z:\'; % Database drive (storage drive)
    system_info.LUT_resolution =  '512f_32w'; %Look-Up Table resolution
    
    signal_info.c = 343; % Speed of sound in metres/sec
    signal_info.Fs = 16000; % Sampling frequency
    signal_info.Nfft = 1024;% Number of fft components
    signal_info.overlap = 0.5;
    signal_info.f_low  = 150;  % Hz
    signal_info.f_high = 8000; % Hz
    signal_info.weight = weight;
    signal_info.L_noise_mask = mask_level; % dB
    signal_info.method = mask_type;
    signal_info.input_filename = [];
    if nargin < 5
        signal_info.L_noise_mask = [];
    end
    system_info.sc = '_'; % Separating character for ascii paths
end

if nargin < 6
    Conv_Type = 'FFT';
end

%% Set values depending on setup type (masker or clean)
if isempty(strfind(lower(signal_info.method), 'masker')) %if not a masker (clean)    
    setup = main_setup;
    levels = -inf;
else %if a masker
    setup = masker_setup;
    levels = signal_info.L_noise_mask;
end

%% Loudspeaker Signals Directory and Recordings Directory
LoudspeakerSignals_Path = cell(length(levels),1);
Recordings_Path = cell(length(levels),1);

for m = 1:length(levels)
    signal_info.L_noise_mask = levels(m);
    LoudspeakerSignals_Path{m} = Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, system_info.LUT_resolution, system_info.Drive );
    Recordings_Path{m} = Results.getRecordingsPath( setup, system_info.LUT_resolution, room_setup, signal_info, system_info.Drive );
end

signal_info.L_noise_mask = levels; % Set values back to originals

%% Load RIR Database file
method = {'new2', 'new', 'old'};
for m = 1:length(method)
    [DB,err] = Room_Acoustics.loadRIRDatabaseFromSetup( main_setup, room_setup, system_info.Drive, method{m});
    if ~err
        break;
    elseif m == length(method)
        error('No matching RIR database file was found. (Generate one and try again)');
    end
end


%% Find Speaker Signals and read to Workspace
fprintf('\n====== Applying Room Impulse Responses to Loudspeaker-Signals ======\n');
fprintf(['            Room Size: ' [strrep(room_setup.Room_Size_txt,'x','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(room_setup.Wall_Absorb_Coeff) '\n']);
fprintf([' Virtual Source Angle: ' num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) '\n']);
fprintf(['       Privacy Method: ' signal_info.method '\n\n']);n=0;
fprintf('\tCompletion: ');

M = length(signal_info.L_noise_mask);
for m = 1:M
    
    files = Tools.getAllFiles( LoudspeakerSignals_Path{m} );
    files = sort(files);
    
    % Continue with only the files that are contained in the original
    % source folder
    files = Tools.keepFilesFromFolder( files, signal_info.speech_filepath);
    if isempty(files), error('No loudspeaker signals found. Have they been generated?'); end    
    
    fileName_prev = '';
    Speaker_Signals = [];
    
    F = length(files);
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
            [fnumflip,fnameflip] = strtok( flip(fileName), SYS.system_info.sc );
            Current_Spkr_Num = str2double( flip( fnumflip ) );
            fileName_curr = flip(sscanf(fnameflip,['%*[' SYS.system_info.sc ']%s']));
            
            if sum(size(y)>1) == 1 && setup.Loudspeaker_Count ~= 1 % If not a multichannel audio file                
                if ~(isempty(fileName_prev) || strcmp( fileName_curr, fileName_prev))
                    Speaker_Signals_ = [];
                end
                
                Speaker_Signals_(Current_Spkr_Num,:) = y;
                if sum( any( Speaker_Signals_, 2 ) ) == SYS.Main_Setup.Loudspeaker_Count
                    Speaker_Signals = Speaker_Signals_;
                end
                
            elseif sum(size(y)>1) >= 1 && setup.Loudspeaker_Count >= 1 % If a multichannel audio file
                Speaker_Signals = y.';
            else
                error(['Failed to read audio file: ' fileName]);
            end
            
            if sum( any( Speaker_Signals, 2 ) ) == setup.Loudspeaker_Count % If we have a complete group of speaker signals
                %% Re-Scale
                % If we have a masker signal then we rescale that signal to
                % its correct level. This is because it was scaled upon
                % saving to audio format to prevent clipping.
                if isfield(signal_info,'clipLevelAdjust')
                    common_scaler = 1/db2mag(signal_info.clipLevelAdjust);
                end
                if ~isempty(strfind(signal_info.method, 'Masker'))
                    if signal_info.L_noise_mask(m) <= 0
                        scaler = (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
                    elseif signal_info.L_noise_mask(m) > 0 % For a positive masker we scaled the signals to save up to a 40db noise masker
                        scaler = (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
                    end
                    Speaker_Signals = Speaker_Signals .* scaler .* common_scaler;
                % The input signal may have also been scaled
                elseif isfield(signal_info,'inputSignalNorm') && signal_info.inputSignalNorm
                    Speaker_Signals = Speaker_Signals .* common_scaler;
                end
                
                %% Parametric Level adjust
                % If we have a parametric array we apply the true
                % ultrasonic reproduction level of that array
                if strcmp( setup.Loudspeaker_Type, 'Parametric' )
                    Speaker_Signals = Speaker_Signals .* setup.Loudspeaker_Object.P1;
                end
                
                %% Perform RIR Convolution with Loudspeaker Signals
                % Bright Zone Signals
                Rec_Sigs_B = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                    Speaker_Signals, DB.RIRs.Bright_RIRs, Conv_Type);
                % Quiet Zone Signals
                Rec_Sigs_Q = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                    Speaker_Signals, DB.RIRs.Quiet_RIRs, Conv_Type);
                
                %% Read original file
                if isempty(strfind(signal_info.method, 'Masker')) % If the signals are not Maskers then read the Original audio
                    try
                        orig = audioread( [LoudspeakerSignals_Path{m}, fileName_curr, system_info.sc, 'Original', fileExt] ); % Assumes same file extension as loudspeaker audio files
                    catch err
                        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                            error('Error reading the Original audio file for the reproduction.'); % Skip unsupported files
                        end
                    end
                else
                    orig=[];
                end
                
                %% Save receiver recordings
                Room_Acoustics.Apply_RIRs.Save_Reverb_Recording( orig, Rec_Sigs_B, Rec_Sigs_Q, signal_info.Fs, Recordings_Path{m}, '', fileName_curr, fileExt );
                
            end
            
            fileName_prev = fileName_curr;
            
        end
        
        [n] = Tools.showTimeToCompletion( ((m-1)*F + f) / (F*M), n);
    end
    
end



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Delete Parallel Pool
delete(para_pool);

end

function setup = swapZones(setup)
Zone = setup.Multizone_Soundfield.Bright_Zone;
setup.Multizone_Soundfield.Bright_Zone = ...
    setup.Multizone_Soundfield.Quiet_Zone;
setup.Multizone_Soundfield.Quiet_Zone = Zone;
end
