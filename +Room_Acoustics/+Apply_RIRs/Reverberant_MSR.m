%% Initialise
clc;
clear;
close all;
fclose all;
tic;
% Start Parallel Pool
para_pool = parpool;
% Kill dropbox
%Tools.Dropbox('kill');
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%% Setup and Path Info
LUT_resolution =  '512f_256w';

Fs = 16000;

%mask_type = 'ZoneWeightMask';
mask_type = 'FlatMask';

%pw_angle = 90;
%pw_angle = 15;
pw_angle = 0;

loudspeakers = 295;
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres
Num_Receivers = 32;

Room_Size = [10 10 10];
%Room_Size = [5 6 4];
%Room_Size = [4 4 3];

Reproduction_Centre = Room_Size ./ 2;%[2.5 3];

Wall_Absorption_Coeff = 1.0;
%Wall_Absorption_Coeff = 0.3;

Drive = 'Z:\';
Spkr_Sig_file_path     = [Drive '+Speaker_Signals\']; % Can be relative or exact
Spkr_Sig_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
Output_file_path_ext = Spkr_Sig_file_path_ext;

%% Load RIR Database file
room = num2str(Room_Size(1));
room = [room 'x' num2str(Room_Size(2)) ];
if length(Room_Size) == 3
    room = [room 'x' num2str(Room_Size(3)) ];
end

room_cent = num2str(Reproduction_Centre(1));
room_cent = [room_cent 'x' num2str(Reproduction_Centre(2)) ];
if length(Room_Size) == 3
    room_cent = [room_cent 'x' num2str(Reproduction_Centre(3)) ];
end

RIR_Database_file_path = [Drive '+Room_Acoustics\+RIR_Database\'];
RIR_Database_file_name = ['RIRs__' num2str(loudspeakers) 'Src_' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab.mat'];;

load([RIR_Database_file_path RIR_Database_file_name]);

ResultsPath = [Drive '+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];

%% Find Speaker Signals and read to Workspace
Input_file_path = [Spkr_Sig_file_path Spkr_Sig_file_path_ext];
files = Tools.getAllFiles(Input_file_path);
files = sort(files);

%Isolate only files of the correct planewave angle and mask type. 
for i=1:length(files)
    ind(i) = (~isempty(findstr(files{i},[num2str(pw_angle) 'pwAngle'])) ...
        && ~isempty(findstr(files{i},['with' mask_type])) ) ...
          || ~isempty(findstr(files{i},'Original'));
end
files = files(ind);

fprintf('\n====== Applying Room Impulse Responses to Loudspeaker-Signals ======\n');
fprintf(['Privacy Weighting: ' mask_type '\n']);n=0;
fprintf('\tCompletion: ');

%parfor_progress( length(files) );

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
    
    
    
    [filePath, fileName, fileExt] = fileparts(files{file});
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
            Rec_Sigs_B = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs(Speaker_Signals, RIRs.Bright_RIRs, 'FFT');
            Rec_Sigs_B = squeeze(sum(Rec_Sigs_B,1));
            % Quiet Zone Signals
            Rec_Sigs_Q = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs(Speaker_Signals, RIRs.Quiet_RIRs, 'FFT');
            Rec_Sigs_Q = squeeze(sum(Rec_Sigs_Q,1));
            
            
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
    %parfor_progress;
end


%parfor_progress(0);

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
%Tools.Dropbox('start');

% Delete Parallel Pool
delete(para_pool);


