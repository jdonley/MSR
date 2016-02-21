function Reverberant_MSR_batchfunc( room_setup, mask_type, loudspeaker_setup, weight, mask_level, Conv_Type )
%% Initialise
tic;
% Start Parallel Pool
para_pool = parpool;
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%% Setup and Path Info
Drive = 'Z:\'; % Database drive (storage drive)

%LUT_resolution =  '512f_256w'; %Look-Up Table resolution
LUT_resolution =  '512f_32w'; %Look-Up Table resolution
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
if nargin < 6
    Conv_Type = 'FFT';
end

sc = '_'; % Separating character for ascii paths

%% Loudspeaker Signals Directory and Recordings Directory
levels = signal_info.L_noise_mask;
LoudspeakerSignals_Path = cell(length(levels),1);
Recordings_Path = cell(length(levels),1);

for m = 1:length(levels)    
    signal_info.L_noise_mask = levels(m);
    LoudspeakerSignals_Path{m} = Broadband_Tools.getLoudspeakerSignalPath( loudspeaker_setup, signal_info, LUT_resolution, Drive );
    Recordings_Path{m} = Results.getRecordingsPath( loudspeaker_setup, LUT_resolution, room_setup, signal_info, Drive );    
end

signal_info.L_noise_mask = levels;

%% Load RIR Database file
% If the signals and setup are for a Masker signal only then reverse the
% zones for the RIR database.
if ~isempty(strfind(signal_info.method, 'Masker'))
    loudspeaker_setup = swapZones(loudspeaker_setup);
end

method = {'new2', 'new', 'old'};
for m = 1:length(method)
    [DB,err] = Room_Acoustics.loadRIRDatabaseFromSetup( loudspeaker_setup, room_setup, Drive, method{m});
    if ~err
        break;
    elseif m == length(method)
        error('No matching RIR database file was found. (Generate one and try again)');
    end
end

% Swap back
if ~isempty(strfind(signal_info.method, 'Masker'))
    loudspeaker_setup = swapZones(loudspeaker_setup);
end

%% Find Speaker Signals and read to Workspace
fprintf('\n====== Applying Room Impulse Responses to Loudspeaker-Signals ======\n');
fprintf(['            Room Size: ' [strrep(room_setup.Room_Size_txt,'x','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(room_setup.Wall_Absorb_Coeff) '\n']);
fprintf([' Virtual Source Angle: ' num2str(loudspeaker_setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) '\n']);
fprintf(['       Privacy Method: ' mask_type '\n\n']);n=0;
fprintf('\tCompletion: ');

M = length(signal_info.L_noise_mask);
for m = 1:M
    
    files = Tools.getAllFiles( LoudspeakerSignals_Path{m} );
    files = sort(files);
    
    fileName_prev = '';
    Speaker_Signals = [];
    
    F = length(files);
    for f = 1:F
        try
            y = audioread(files{f});
        catch err
            if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                continue; % Skip unsupported files
            end
        end
        
        
        
        [~, fileName, fileExt] = fileparts(files{f});
        if isempty(strfind(fileName,'Original')) % Make sure the file being read isn't an original file
            
            % Get the file number and file name
            [fnumflip,fnameflip] = strtok( flip(fileName), sc );
            Current_Spkr_Num = str2double( flip( fnumflip ) );
            fileName_curr = flip(sscanf(fnameflip,['%*[' sc ']%s']));
            
            if ~(isempty(fileName_prev) || strcmp( fileName_curr, fileName_prev))
                Speaker_Signals = [];
            end
            
            Speaker_Signals(Current_Spkr_Num,:) = y;
            
            if sum( any( Speaker_Signals, 2 ) ) == loudspeaker_setup.Loudspeaker_Count % If we have a complete group of speaker signals
                %% Re-Scale
                % If we have a masker signal then we rescale that signal to
                % its correct level. This is because it was scaled upon
                % saving to audio format to prevent clipping.
                if ~isempty(strfind(signal_info.method, 'Masker'))
                    if signal_info.L_noise_mask(m) <= 0
                        scaler = (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
                    elseif signal_info.L_noise_mask(m) > 0 % For a positive masker we scaled the signals to save up to a 40db noise masker
                        scaler = (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
                    end
                    Speaker_Signals = Speaker_Signals .* scaler;
                end
                
                %% Parametric Level adjust
                % If we have a parametric array we apply the true
                % ultrasonic reproduction level of that array
                if strcmp( loudspeaker_setup.Loudspeaker_Type, 'Parametric' )
                   Speaker_Signals = Speaker_Signals .* loudspeaker_setup.Loudspeaker_Object.P1;
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
                        orig = audioread( [LoudspeakerSignals_Path{m}, fileName_curr, sc, 'Original', fileExt] ); % Assumes same file extension as loudspeaker audio files
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
