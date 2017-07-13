function Calibrate_MSR_Loudspeaker_Signals(SYS)
clear;

%% Setup Information
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end

%Flip loudspeaker order (effectively flips entire setup) (false if not needed)
if isfield(SYS.system_info,'MirrorSetup'), MirrorSetup = SYS.system_info.MirrorSetup;
else MirrorSetup = false; end

% If a realworld recording is not specified in the system then abort
if ~any(strcmpi(strrep(SYS.signal_info.recording_type,'-',''),'realworld')), delete(gcp('nocreate')); return; end

%%
% for ss = 1:numel(SYS.Main_Setup)
%     SYS_ = SYS;
%     SYS_.Main_Setup = SYS_.Main_Setup(ss);

N = length(SYS.signal_info.methods_list_clean(SYS.signal_info.methods_list_clean>=1)) ...
    + length(SYS.signal_info.methods_list_masker(SYS.signal_info.methods_list_masker>=1));
paired = isfield(SYS.signal_info,'methods_list_paired') && SYS.signal_info.methods_list_paired;

for typ = 1:N
    
    SYS.signal_info.method = SYS.signal_info.methods_list{typ};
    SYS_ = SYS;
    
    tmp_levels = SYS_.signal_info.L_noise_mask;
    
    if any( typ == SYS_.signal_info.methods_list_clean )
        Setup = SYS_.Main_Setup(typ == SYS_.signal_info.methods_list_clean);
        if paired
            SYS_.Masker_Setup = SYS_.Masker_Setup(typ == SYS_.signal_info.methods_list_clean);
        end
        noise_levels = -inf;
    elseif any( typ == SYS_.signal_info.methods_list_masker )
        Setup = SYS_.Masker_Setup(typ == SYS_.signal_info.methods_list_masker);
        if paired
            SYS_.Main_Setup = SYS_.Main_Setup(typ == SYS_.signal_info.methods_list_masker);
        end
        noise_levels = SYS_.signal_info.L_noise_mask;
    else
        error(['The method ''' SYS_.signal_info.method ''' is not set as a clean or masker signal'])
    end
    
    
    fprintf('\n===== Calibrating Loudspeaker Signals =====\n');
    fprintf('\tSignal Type: %s\n', SYS_.signal_info.method);
    fprintf('\t Completion: '); n=0; tic;
    M = length(noise_levels);
    
    %% Load Fiters
    filter_location = Tools.getAllFiles([SYS_.system_info.Drive SYS_.system_info.FilterData_dir]);
    filter_location(~contains( filter_location, 'EQ'))=[];
    filter_location = sort(filter_location);
    if isempty(filter_location)
       Tools.simpleWarning('No calibration filters were found. Skipping calibration.');
       return;
    end
    filts = load( filter_location{1} ); % 1st element should be the newest (most recent) set of filters
    if MirrorSetup
        filts.EQ = flip(filts.EQ,2);
    end
    
    for m = 1:M
        noise_mask = noise_levels(m);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SYS_.signal_info.L_noise_mask = noise_mask; % dB
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% Get Paths and files
        [Spkr_path,~,~,~,~,~,path_ext] = ...
            Broadband_Tools.getLoudspeakerSignalPath( Setup, SYS_.signal_info, SYS_.system_info.LUT_resolution, SYS_.system_info.Drive, 'new');
        
        files = Tools.getAllFiles(Spkr_path);
        files = sort(files);
        
        isMasker = contains(lower(SYS_.signal_info.method), 'masker'); %if a masker
        % Continue with only the files that are contained in the original
        % source folder unless a Masker signal is being processed (these are
        % named as Maskers)
        if ~isMasker
            files = Tools.keepFilesFromFolder( files, SYS_.signal_info.speech_filepath);
        end
        if isempty(files), error('No loudspeaker signals found. Have they been generated?'); end
        
        %%
        fileName_prev = '';
        Speaker_Signals = [];
        Speaker_Signals_ = [];
        
        F=length(files);
        for f = 1:F
            
            [~, fileName, fileExt] = fileparts(files{f});
            if ~contains(  fileName, 'Original'  ) % If not an Original audio file
                try
                    y = audioread(files{f});
                catch err
                    if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                        continue; % Skip unsupported files
                    end
                end
                
                % Get the file number and file name
                [fnumflip,fnameflip] = strtok( flip(fileName), SYS_.system_info.sc );
                Current_Spkr_Num = str2double( flip( fnumflip ) );
                fileName_curr = flip(sscanf(fnameflip,['%*[' SYS_.system_info.sc ']%s']));
                
                if sum(size(y)>1) == 1 % If not a multichannel audio file
                    if ~(isempty(fileName_prev) || strcmp( fileName_curr, fileName_prev))
                        Speaker_Signals_ = [];
                    end
                    
                    Speaker_Signals_(Current_Spkr_Num,:) = y;
                    if sum( any( Speaker_Signals_, 2 ) ) == Setup.Loudspeaker_Count
                        Speaker_Signals = Speaker_Signals_;
                    end
                    
                elseif sum(size(y)>1) > 1 % If a multichannel audio file
                    Speaker_Signals = y.';
                else
                    error(['Failed to read audio file: ' fileName]);
                end
                
                if sum( any( Speaker_Signals, 2 ) ) == Setup.Loudspeaker_Count % If we have a complete group of speaker signals
                    Speaker_Signals = Speaker_Signals.';
                    Speaker_Signals = flip(Speaker_Signals,2); % Speakers are in reverse order
                    
                    %% Read original file
                    if isempty(strfind(SYS_.signal_info.method, 'Masker')) % If the signals are not Maskers then read the Original audio
                        try
                            orig = audioread( [Spkr_path, fileName_curr, SYS_.system_info.sc, 'Original', fileExt] ); % Assumes same file extension as loudspeaker audio files
                        catch err
                            if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                                error('Error reading the Original audio file for the reproduction.'); % Skip unsupported files
                            end
                        end
                    else
                        orig=[];
                    end
                    
                    %% Upsample to the filter sampling frequency
                    up_factor = filts.fs/SYS_.signal_info.Fs;
                    [b,a] = butter( 9, 1/up_factor );
                    loudspeaker_signals_upsampled = zeros(size(Speaker_Signals,1)*up_factor, Setup.Loudspeaker_Count);
                    for spkr = 1:Setup.Loudspeaker_Count
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
                    
                    output_dir = [SYS_.system_info.Drive SYS_.system_info.Calibrated_Signals_dir path_ext];
                    
                    if ~exist(output_dir,'dir'); mkdir(output_dir); end
                    
                    if sum(size(y)>1) == 1 % If not a multichannel audio file
                        for spkr = 1:Setup.Loudspeaker_Count
                            audiowrite([output_dir fileName_curr SYS_.system_info.sc 'Upsampled' SYS_.system_info.sc num2str(spkr) fileExt], ...
                                loudspeaker_signals_calibrated(:,spkr), filts.fs);
                        end
                    elseif sum(size(y)>1) > 1 % If a multichannel audio file
                        audiowrite([output_dir fileName_curr SYS_.system_info.sc 'Upsampled' fileExt], ...
                            loudspeaker_signals_calibrated, filts.fs);
                    end
                    save([output_dir fileName_curr SYS_.system_info.sc 'Upsampled' '.mat'], 'scale_value');
                    
                    if ~isempty(orig)
                        audiowrite([output_dir fileName_curr SYS_.system_info.sc 'Original' fileExt], orig_calibrated, filts.fs );
                    end
                    
                end
                
                fileName_prev = fileName_curr;
            end
            n = Tools.showTimeToCompletion( (f + (m-1)*F) / (F*M), n);
        end
    end
    
    fprintf('\nCompleted\n\n');
    
    %Return values to normal
    SYS_.signal_info.L_noise_mask = tmp_levels;
    %     end
end




