function Reverberant_MSR_Analysis_batchfunc( SYS_or_setups, room_setup, signal_types, weight, mask_level, pesqNumber, RecordingType)

SYS_type = 'Current_Systems.SR_System';

%% Initialise
tic;
% Start Parallel Pool
para_pool = parpool;
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%%
if nargin == 1
    if isa(SYS_or_setups,SYS_type)
        SYS = SYS_or_setups;
        Mains = num2cell(SYS.Main_Setup);
        Maskers = num2cell(SYS.Masker_Setup);
        setups = {Mains{:},Maskers{:}};
        room_setup = SYS.Room_Setup;
        signal_types = SYS.signal_info.methods_list;
        signal_info = SYS.signal_info;
        system_info = SYS.system_info;
        analysis_info = SYS.analysis_info;
        isHybrid = size(SYS.signal_info.methods_list_masker,1)>1;
    else
        error(['Single input argument must be of type: ' SYS_type]);
    end
elseif nargin > 1
    setups = SYS_or_setups;
    system_info.Drive = 'Z:\'; % Database drive (storage drive)
    system_info.LUT_resolution =  '512f_32w'; %Look-Up Table resolution
    signal_info.c = 343; % Speed of sound in metres/sec
    signal_info.Fs = 16000; % Sampling frequency
    signal_info.Nfft = 1024;% Number of fft components
    signal_info.overlap = 0.5;
    signal_info.f_low  = 150;  % Hz
    signal_info.f_high = 8000; % Hz
    signal_info.f_low_meas = 100; % Hz %Minimum loudspeaker response
    signal_info.f_high_meas = 7000; % Hz %Maximum frequency with accurate response at given sampling rate
    signal_info.weight = weight;
    signal_info.L_noise_mask = mask_level; % dB
    signal_info.input_filename = [];
    analysis_info.Measures = {'PESQ','STOI','SNR'};
    isHybrid = length(signal_types) >= 3;
    if nargin < 7
        RecordingType = 'simulated';
    end
    signal_info.recording_type = RecordingType;
    if nargin < 5
        signal_info.L_noise_mask = [];
    end
    system_info.sc = '_'; % Separating character for ascii paths
end

if nargin < 6
    pesqNumber = 0;
end
Measures = analysis_info.Measures;

signal_type = [signal_types{2:end}];


signal_info.recording_type = strrep(signal_info.recording_type,'-','');
isSimulated = strcmpi(signal_info.recording_type,'simulated');
isRealworld = strcmpi(signal_info.recording_type,'realworld');
isSweep     = isempty(signal_type) && any(strcmp(analysis_info.Measures,'SPL'));


if isHybrid
    signal_type = ['Hybrid' signal_type];
end
if isSweep
    signal_type = 'Sweep'; % An exponential sine sweep is what is used to measure SPL with this framework
end
signal_info.method = signal_type;


%% Obtain Recordings and Results Directory Path
ResultsPath = Results.getResultsPath( setups{1}, system_info.LUT_resolution, room_setup, signal_info, system_info.Drive );
if ~exist(ResultsPath,'dir'); mkdir(ResultsPath); end
save([ResultsPath 'Setups.mat'], 'setups');

levels = signal_info.L_noise_mask;
weight = signal_info.weight;
Recordings_Path = cell(length(levels),length(setups));

if isSimulated
    s_1 = 1;
elseif isRealworld
    s_1 = 2;
    if isSweep
        s_1 = 1;
    end 
end

for s = s_1:length(setups)
    signal_info.method = signal_types{s};
    if isempty(strfind( signal_types{s}, 'Parametric' )) % If not parametric and more than one loudspeaker
        signal_info.weight = weight;
    else % If parametric and/or only one loudspeaker
        signal_info.weight = 1;
    end
    for m = 1:length(levels)
        if any(s == signal_info.methods_list_clean)
            signal_info.L_noise_mask = -Inf; %Masker level is non-existent for clean setup (speech setup)
        elseif any(s == signal_info.methods_list_masker)
            signal_info.L_noise_mask = levels(m);
        else
            error('Setup index not defined in ''methods_list_clean'' or ''methods_list_masker''.')
        end
        if isSimulated
            Recordings_Path{m,s} = Results.getRecordingsPath( setups{s}, system_info.LUT_resolution, room_setup, signal_info, system_info.Drive );
        elseif isRealworld
            Recordings_Path{m,s} = Hardware_Control.getRealRecordingsPath( setups{1}, system_info.LUT_resolution, room_setup, signal_info, system_info.Drive );
        end
    end
end
signal_info.L_noise_mask = levels;
signal_info.weight = weight;
signal_info.method = signal_type;


%% Delete previous results file
Results.deleteResultsFile( ResultsPath, analysis_info.Measures);

%% Start Evaluation Loop
if isSimulated
    EnvType = 'Simulated';
elseif isRealworld
    EnvType = 'Real-World';
end
fprintf(['\n====== Analysing ' EnvType ' Reverberant Signals ======\n']);
fprintf(['            Room Size: ' [strrep(sprintf(strrep(repmat('%g',1,length(room_setup.Room_Size)),'g%','g %'),room_setup.Room_Size),' ','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(room_setup.Wall_Absorb_Coeff) '\n']);
fprintf([' Virtual Source Angle: ' num2str(setups{1}.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) '\n']);
fprintf(['    Privacy Weighting: ' signal_info.method '\n\n']);n=0;
fprintf('\t Completion: ');


%% Find Speaker Signals and read to Workspace
M = length(signal_info.L_noise_mask);
S = length(setups);
for m = 1:M
    
    files = {[],[]};
    for s=s_1:S
        files_ = Tools.getAllFiles( Recordings_Path{m,s} );
        % Continue with only the files that are contained in the original
        % source folder
        if s==s_1
            files_ = Tools.keepFilesFromFolder( files_, signal_info.speech_filepath);
        end
        % Warn if there are no recordings to analyse and then return from the function
        if isempty(files_)
            wrnCol = [255,100,0]/255;
            cprintf(wrnCol, 'There was a problem reading associated recording files.\n');
            cprintf(wrnCol, ['Do ''' signal_info.recording_type ''' recordings exist?\n']);
            cprintf(wrnCol, ['Skipping ''' signal_info.recording_type ''' analysis procedure.\n']);
            delete(gcp('nocreate')); return;
        end
        files{:,s} = sort(files_);
    end
    
    fileName_prev = '';
    Rec_Bright = [];
    Rec_Quiet = [];
    Fs = signal_info.Fs;
    
    % Load maskers
    for s=s_1:S
        if any(s == signal_info.methods_list_masker) || ...
            (all(signal_info.methods_list_masker==0) && any(strcmp(analysis_info.Measures,'SPL')))
            fileName={};
            i=(all(signal_info.methods_list_masker==0) && any(strcmp(analysis_info.Measures,'SPL')))*1;
            fpos = [1 (s_1+1+i)];
            for f=fpos
                [~, fileName{s,f}, ~] = fileparts(files{s}{f});
            end
            % These should be maskers unless realworld recording
            Rec_Bright_{s} = load(files{s}{fpos(1)});
            Rec_Quiet_{s} = load(files{s}{fpos(2)});
            if isfield(Rec_Bright_{s},'fs')
                Fs = Rec_Bright_{s}.fs;
            end
        end
    end
    %     if S == 3 % If two maskers are found
    %         [Rec_Bright_, Rec_Quiet_] = adjustHybridMaskers( ...
    %             Rec_Bright_, Rec_Quiet_, setups, signal_info);
    %     end
    
    F = size(files{s_1},1);
    for file = 1:F
        
        [~, fileName{s_1,file}, ~] = fileparts(files{s_1}{file});
        
        
        if isempty(strfind(fileName{s_1,file},'Original')) % Make sure the file being read isn't an original file
            
            % Get the file number and file name
            [Ztypeflip,Stypeflip] = strtok( flip(fileName{s_1,file}), system_info.sc );
            SignalName = flip( Stypeflip );
            ZoneType = flip( Ztypeflip );
            
            
            if ~(isempty(fileName_prev) || strcmp( SignalName, fileName_prev))
                Rec_Bright = [];
                Rec_Quiet = [];
            end
            
            if strcmp('Bright',ZoneType)
                Rec_Bright = [];
                for s = signal_info.methods_list_clean.'
                    if ~isempty(files{s})
                        Rec_Bright_{s} = load(files{s}{file});
                    else
                        Rec_Bright_{s_1} = load(files{s_1}{file}); break;
                    end
                end
                if s_1==2, Rec_Bright_{s_1}.Rec_Sigs_B = Rec_Bright_{s_1}.Rec_Sigs_B'; end %TODO: Fix the recording so the dimensions are in the correct place.
                sLB = size( Rec_Bright_{s_1}.Rec_Sigs_B,2); %signal Length Bright
                for s=s_1:S
                    Rec_Bright(:,:,s) = Rec_Bright_{s}.Rec_Sigs_B(:,1:sLB);
                end
                Rec_Bright = sum( Rec_Bright, 3 );
                [~,Iforce]=sort(size(Rec_Bright)); 
                Rec_Bright=permute(Rec_Bright,Iforce);% Force smaller dimension first
            elseif strcmp('Quiet',ZoneType)
                Rec_Quiet = [];
                for s = signal_info.methods_list_clean.'
                    if ~isempty(files{s})
                        Rec_Quiet_{s} = load(files{s}{file});
                    else
                        Rec_Quiet_{s_1} = load(files{s_1}{file}); break;
                    end
                end
                if s_1==2, Rec_Quiet_{s_1}.Rec_Sigs_Q = Rec_Quiet_{s_1}.Rec_Sigs_Q'; end; %TODO: Fix the recording so the dimensions are in the correct place.
                sLQ = size( Rec_Quiet_{s_1}.Rec_Sigs_Q,2); %signal Length Quiet
                for s=s_1:S
                    Rec_Quiet(:,:,s) = Rec_Quiet_{s}.Rec_Sigs_Q(:,1:sLQ);
                end
                Rec_Quiet = sum( Rec_Quiet, 3 );
                [~,Iforce]=sort(size(Rec_Quiet)); 
                Rec_Quiet=permute(Rec_Quiet,Iforce);% Force smaller dimension first
            else
                error('Error reading the recordings. The type of zone read from the file name is not supported');
            end
            
            if all( [any(Rec_Bright(:)), any(Rec_Quiet(:))] ) % If we have two complete group of receiver signals for each zone
                
                % Read original file
                for file_orig = 1:F
                    if ~isempty( strfind(files{s_1}{file_orig}, [SignalName 'Original']) )
                        break
                    end
                end
                try
                    orig = audioread( files{s_1}{file_orig} );
                catch err
                    if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                        continue; % Skip unsupported files
                    end
                end
                
                
                % BEGIN Downsample realworld recordings
                if Fs ~= signal_info.Fs
                    down_rate = Fs / signal_info.Fs ;
                    Rec_Bright_down = zeros(ceil(size(Rec_Bright).*[1 1/down_rate]));
                    Rec_Quiet_down = zeros(ceil(size(Rec_Quiet).*[1 1/down_rate]));
                    for  r = 1:size(Rec_Bright,1)
                        Rec_Bright_down(r,:) = ...
                            decimate( Rec_Bright(r,:), down_rate );
                        Rec_Quiet_down(r,:) = ...
                            decimate( Rec_Quiet(r,:), down_rate );
                    end
                    Rec_Bright = Rec_Bright_down;
                    Rec_Quiet  = Rec_Quiet_down ;
                    if ~isempty(orig)
                        orig = decimate( orig, down_rate );
                    end
                end
                % END Downsample realworld recordings
                
                if ~any(cell2mat(strfind(upper(Measures),'SUPPRESSION')))
                    % BEGIN filter signals to acceptable measurement frequency range
                    [b,a] = butter(6, [signal_info.f_low_meas signal_info.f_high_meas] ./ (signal_info.Fs/2) );
                    %[b,a] = cheby1(6 [signal_info.f_low_meas signal_info.f_high_meas] ./ (signal_info.Fs/2) );
                    orig = filter(b,a,orig);
                    Rec_Bright = filter(b,a,Rec_Bright')';
                    Rec_Quiet  = filter(b,a,Rec_Quiet')';
                    % END filter signals
                    
                    % BEGIN power normalise to bright zone level
                    adjScale = [];
                    for r = 1:size(Rec_Bright,1)
                        [~,adjScale(r)] = Broadband_Tools.power_norm(orig, Rec_Bright(r,:), signal_info.Fs, [signal_info.f_low_meas signal_info.f_high_meas]);
                    end
                    Rec_Bright = Rec_Bright * mean(adjScale);
                    Rec_Quiet  = Rec_Quiet * mean(adjScale);
                    % END power normalise
                    
                    % BEGIN Resize the original speech signal and
                    % Align it with the reverberant signals.
                    orig(length(orig):size(Rec_Bright,2))=0; % Resize the original signal because the reverberant signal will be longer
                    if (length(orig) ~= length(Rec_Bright)) || (length(orig) ~= length(Rec_Quiet))
                        error('Size of the original signal does not match the reproduced signal!');
                    end
                    
                    %c_speed = 343;%343m/s speed of sound in air
                    %max_delay = speaker_radius*2 / c_speed * signal_info.Fs;
                    max_delay = signal_info.Fs / 2;
                    Original = zeros( size(Rec_Bright,1), 2, ...
                        max( [length(orig), size(Rec_Bright,2), size(Rec_Quiet,2)] ) + max_delay );
                    Rec_BrightTMP = zeros(size(Rec_Bright,1),size(Original,3));
                    Rec_QuietTMP  = zeros(size(Rec_Quiet ,1),size(Original,3));
                    for r = 1:size(Rec_Bright,1)
                        [tmp1, tmp2] = alignsignals( orig, Rec_Bright(r,:), max_delay);
                        [tmp3, tmp4]  = alignsignals( orig, Rec_Quiet(r,:),  max_delay);
                        Original(r,1,1:numel(tmp1)) = tmp1;
                        Rec_BrightTMP(r,1:numel(tmp2)) = tmp2;
                        Original(r,2,1:numel(tmp3)) = tmp3;
                        Rec_QuietTMP(r,1:numel(tmp4))  = tmp4;
                    end
                    Rec_Bright = Rec_BrightTMP;
                    Rec_Quiet = Rec_QuietTMP;
                    % END resize and align
                end
                
                % BEGIN Calculate and save results
                %%% ANALYSIS
                % Speech Quality
                if any(cell2mat(strfind(upper(Measures),'PESQ')))
                    % PESQ - Perceptual Evaluation of Speech Quality
                    Room_Acoustics.Apply_RIRs.Save_Reverb_PESQ_Result( Original, Rec_Bright, signal_info.Fs, signal_info.L_noise_mask(m), ResultsPath, [], SignalName, pesqNumber );
                end
                
                % Speech Intelligibility
                if any(cell2mat(strfind(upper(Measures),'STOI')))
                    % STOI - Short-Time Objective Intelligibility
                    Room_Acoustics.Apply_RIRs.Save_Reverb_STOI_Result( Original, Rec_Bright, Rec_Quiet, signal_info.Fs, signal_info.L_noise_mask(m), ResultsPath, [], SignalName );
                end
                if any(cell2mat(strfind(upper(Measures),'STI')))
                    % STI - Speech Transmission Index
                    Room_Acoustics.Apply_RIRs.Save_Reverb_STI_Result( Original, Rec_Bright, Rec_Quiet, signal_info.Fs, ResultsPath, [], SignalName );
                end
                
                if any(cell2mat(strfind(upper(Measures),'SNR')))
                    % SNR - Signal to Noise Ratio
                    Room_Acoustics.Apply_RIRs.Save_Reverb_SNR_Result( Original, Rec_Bright, Rec_Quiet, signal_info.Fs, signal_info.L_noise_mask(m), ResultsPath, [], SignalName );
                end
                
                if any(cell2mat(strfind(upper(Measures),'SUPPRESSION')))
                    % Interference Suppression
                    Room_Acoustics.Apply_RIRs.Save_Reverb_Suppression_Result( Rec_Bright_{end}.Rec_Sigs_B(:,1:sLB), Rec_Bright, signal_info.Fs, ...
                        analysis_info.Nfft, [analysis_info.f_low, analysis_info.f_high], ResultsPath, [], SignalName );
                end
                
                if any(cell2mat(strfind(upper(Measures),'SPL')))
                    % Sound Pressure Levels
                    Room_Acoustics.Apply_RIRs.Save_Reverb_SPLs_Result( Rec_Bright, Rec_Quiet, signal_info.Fs, ...
                        analysis_info.Nfft, [analysis_info.f_low, analysis_info.f_high], ResultsPath, [], SignalName, SYS );
                end
                % END calc and save results
                
            end
            
            fileName_prev = SignalName;
            
        end
        
        n = Tools.showTimeToCompletion( ((m-1)*F + file)/ (F*M), n);
        
    end
end

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Delete Parallel Pool
delete(para_pool);

end

function [Rec_B, Rec_Q] = adjustHybridMaskers(Rec_B, Rec_Q, Setups, SigInfo) % Assumes two maskers only
% Find cutoff from multizone setup
S = Setups{1};
BZ = S.Multizone_Soundfield.Bright_Zone;
QZ = S.Multizone_Soundfield.Quiet_Zone;
R_ = max( [BZ.Radius_q + BZ.Origin_q.Distance; ...
    QZ.Radius_q + QZ.Origin_q.Distance;  ]);
phiL_rad = S.Speaker_Arc_Angle / 180 * pi;
f_cutoff = SigInfo.c * (S.Loudspeaker_Count - 1) / (2 * R_ * phiL_rad);
f_cutoff = mean(f_cutoff .* [1, 2]);

% Match two signals
% Adjust second masker to suit first
[newsig, adjVal] = Broadband_Tools.power_norm( Rec_Q{2}.Rec_Sigs_Q', Rec_Q{3}.Rec_Sigs_Q', SigInfo.Fs, f_cutoff );
Rec_Q{3}.Rec_Sigs_Q = newsig';
Rec_B{3}.Rec_Sigs_B = Rec_B{3}.Rec_Sigs_B * adjVal; % Same adjustment to bright second masker

end

