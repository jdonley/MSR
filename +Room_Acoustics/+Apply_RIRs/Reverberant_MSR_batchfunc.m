function Reverberant_MSR_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, speech_setup, mask_level, Conv_Type )
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
Fs = 16000; %Sampling Frequency

pw_angle = speech_setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle;
        
Num_Receivers = 32; %Number of recording points (microphones)

Reproduction_Centre = Room_Size ./ 2;

if nargin < 5
   mask_level = [];
end
if nargin < 6
    Conv_Type = 'FFT';
end

%%
Spkr_Sig_file_path     = [Drive '+Speaker_Signals\']; % Can be relative or exact
Spkr_Sig_file_path_ext = ['+' num2str(speech_setup.Radius*2) 'm_SpkrDia\+' num2str(speech_setup.Loudspeaker_Count) 'Spkrs_' num2str(speech_setup.Speaker_Arc_Angle) 'DegArc_LUT_' LUT_resolution '\'];
Output_file_path_ext = Spkr_Sig_file_path_ext;

%% Load RIR Database file
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

RIR_Database_file_path = [Drive '+Room_Acoustics\+RIR_Database\'];
RIR_Database_file_name = ['RIRs__' num2str(speech_setup.Loudspeaker_Count) 'Src_' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab.mat'];

load([RIR_Database_file_path RIR_Database_file_name]);

ResultsPath = [Drive '+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];

%% Find Speaker Signals and read to Workspace
Input_file_path = [Spkr_Sig_file_path Spkr_Sig_file_path_ext];
files = Tools.getAllFiles(Input_file_path);

%Isolate only files of the correct planewave angle, mask type and mask level.
for i=1:length(files)
    tmp = 1;
    for m = 1:length(mask_level)
        tmp(m) = ~isempty(strfind(files{i},['_' num2str(mask_level(m)) 'dB']));
    end
    tmp = any(tmp);
    ind(i) =   (~isempty(strfind(files{i},['_' num2str(pw_angle) 'pwAngle'])) ...
                && ~isempty(strfind(files{i},['with' mask_type])) ...
                && tmp ) ...
            || (~isempty(strfind(files{i},['_' num2str(pw_angle) 'pwAngle'])) ...
                && isempty(strfind(files{i},'with')) ...
                && strcmp(mask_type,'NoMask')) ...
            || ~isempty(strfind(files{i},'Original'));
end
files = files(ind);
files = sort(files);

fprintf('\n====== Applying Room Impulse Responses to Loudspeaker-Signals ======\n');
fprintf(['            Room Size: ' [strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(Wall_Absorption_Coeff) '\n']);
fprintf(['      Planewave Angle: ' num2str(pw_angle) '\n']);
fprintf(['    Privacy Weighting: ' mask_type '\n\n']);n=0;
fprintf('\tCompletion: ');


fileName_prev = '';
Speaker_Signals = [];
for file = 1:length(files)
    try
        y = audioread(files{file});
    catch err
        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
            continue; % Skip unsupported files
        end
    end
    
    
    
    [~, fileName, fileExt] = fileparts(files{file});
    if isempty(strfind(fileName,'Original')) % Make sure the file being read isn't an original file
        
        Num_ = strfind(fileName,'__')+2;
        Num_p = Num_(3);
        
        if ~(strcmp(fileName_prev, '') || strcmp(strrep(fileName(1:Num_p-1),'_',''), fileName_prev))
            Speaker_Signals = [];
        end
        
        Current_Spkr_Num = sscanf(fileName(Num_p-1:end),'%*[^0-9]%d');
        Speaker_Signals(Current_Spkr_Num,:) = y;
        
        if all( any( Speaker_Signals, 2 ) ) % If we have a complete group of speaker signals
            
            % Perform RIR Convolution with Loudspeaker Signals
            % Bright Zone Signals
            Rec_Sigs_B = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                                Speaker_Signals, RIRs.Bright_RIRs, Conv_Type);
            % Quiet Zone Signals
            Rec_Sigs_Q = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                                Speaker_Signals, RIRs.Quiet_RIRs, Conv_Type);
            
            
            % Read original file
            for file_orig = 1:length(files)
                if ~isempty( strfind(files{file_orig}, [fileName(1:Num_(1)-3) '__Original_']) )
                    break
                end
            end
            try
                orig = audioread( files{file_orig} );
            catch err
                if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                    continue; % Skip unsupported files
                end
            end
             
            % Save receiver recordings
            Room_Acoustics.Apply_RIRs.Save_Reverb_Recording( orig, Rec_Sigs_B, Rec_Sigs_Q, Fs, ResultsPath, Output_file_path_ext, fileName, fileExt );
                  
        end
        
        fileName_prev = strrep(fileName(1:Num_p-1),'_','');
    
    end
    
    n = Tools.showTimeToCompletion( file/length(files), n);
end



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Delete Parallel Pool
delete(para_pool);


