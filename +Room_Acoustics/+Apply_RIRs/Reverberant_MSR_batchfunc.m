function Reverberant_MSR_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, pw_angle, mask_level)
%% Initialise
%clc;
%clear;
%close all;
%fclose all;
tic;
% Start Parallel Pool
para_pool = parpool;
% Kill dropbox
%Tools.Dropbox('kill');
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%% Setup and Path Info
LUT_resolution =  '512f_256w'; %Look-Up Table resolution
Fs = 16000; %Sampling Frequency
loudspeakers = 295; %Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres
Time_Delay = false;
Num_Receivers = 32; %Number of recording points (microphones)
if nargin < 5
   mask_level = [];
end

%%
% % ROOM 1
% % Anechoic
%Room_Size = [10 10 10]; %Anechoic
%Wall_Absorption_Coeff = 1.0;

% % ROOM 2
% % Small Open Plan Office
% %Room_Size = [5 6 4];  %Small Open Plan Office
% %Wall_Absorption_Coeff = 0.3;

% % ROOM 3
% % Medium Open Plan Office
% Room_Size = [8 10 3];   %Medium Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 4
% % Cafe / Restaurant
%Room_Size = [9 14 3];   %Cafe/Restaurant
%Wall_Absorption_Coeff = 0.3;

% % ROOM 5
% % Small Open Plan Office
%Room_Size = [4 9 3];   %Small Open Plan Office
%Wall_Absorption_Coeff = 0.3;

%%
% % Setup and Privacy Scheme 1
%mask_type = 'FlatMask';
%pw_angle = 0;

% % Setup and Privacy Scheme 2
%mask_type = 'FlatMask';
%pw_angle = 15;

    % % Setup and Privacy Scheme 3
    %mask_type = 'ZoneWeightMask';
    %pw_angle = 15;

% % Setup and Privacy Scheme 4
%mask_type = 'FlatMask';
%pw_angle = 90;

    % % Setup and Privacy Scheme 5
    %mask_type = 'ZoneWeightMask';
    %pw_angle = 90;

%%
Reproduction_Centre = Room_Size ./ 2;

Drive = 'Z:\';
Spkr_Sig_file_path     = [Drive '+Speaker_Signals\']; % Can be relative or exact
Spkr_Sig_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
Output_file_path_ext = Spkr_Sig_file_path_ext;

%% Load RIR Database file
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

RIR_Database_file_path = [Drive '+Room_Acoustics\+RIR_Database\'];
RIR_Database_file_name = ['RIRs__' num2str(loudspeakers) 'Src_' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab.mat'];;

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

        
        ind(i) = ~isempty(strfind(files{i},'sinusoid'));


end
files = files(ind);
files = sort(files);

fprintf('\n====== Applying Room Impulse Responses to Loudspeaker-Signals ======\n');
fprintf(['            Room Size: ' [strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(Wall_Absorption_Coeff) '\n']);
fprintf(['      Planewave Angle: ' num2str(pw_angle) '\n']);
fprintf(['    Privacy Weighting: ' mask_type '\n\n']);n=0;
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
            
            %Time delay (group delay) loudspeaker signals for better
            %reproduction in target brigh zone
            if Time_Delay
                load('+Room_Acoustics\Time_Delay_TFs.mat');
                TD_TF=[]; TD_TF(1,:,:) = Time_Delay_TFs';
                Speaker_Signals = squeeze( ...
                    Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                    Speaker_Signals, TD_TF, 'FFT'));
            end
            
            % Perform RIR Convolution with Loudspeaker Signals
            % Bright Zone Signals
            Rec_Sigs_B = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                                Speaker_Signals, RIRs.Bright_RIRs, 'FFT');
            Rec_Sigs_B = squeeze(sum(Rec_Sigs_B,1));
            % Quiet Zone Signals
            Rec_Sigs_Q = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs( ...
                                Speaker_Signals, RIRs.Quiet_RIRs, 'FFT');
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


