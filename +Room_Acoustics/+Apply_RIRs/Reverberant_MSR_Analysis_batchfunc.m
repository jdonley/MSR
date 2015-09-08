function Reverberant_MSR_Analysis_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, pw_angle, mask_level)
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

Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];

%% Load RIR Database file
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

Drive = 'Z:\';
ResultsPath = ['+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];

%% Find Speaker Signals and read to Workspace
Input_file_path = [Drive ResultsPath Output_file_path_ext];
files = Tools.getAllFiles(Input_file_path);
files = sort(files);

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

fprintf('\n====== Analysing Simulated Reverberant Signals ======\n');
fprintf(['            Room Size: ' [strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(Wall_Absorption_Coeff) '\n']);
fprintf(['      Planewave Angle: ' num2str(pw_angle) '\n']);
fprintf(['    Privacy Weighting: ' mask_type '\n\n']);n=0;
fprintf('\tCompletion: ');

%parfor_progress( length(files) );

fileName_prev = '';
Rec_Bright = [];
Rec_Quiet = [];
for file = 1:length(files)
   
    [filePath, fileName, fileExt] = fileparts(files{file});
    
    if isempty(strfind(fileName,'Original')) % Make sure the file being read isn't an original file
        
        Num_ = strfind(fileName,'__')+2;
        Num_p = Num_(3);
        
        if ~(strcmp(fileName_prev, '') || strcmp(strrep(fileName(1:Num_p-1),'_',''), fileName_prev))
            Rec_Bright = [];
            Rec_Quiet = [];
        end
        
        Zone_Type = sscanf([' ' fileName(Num_p:end)],'%*[^_]_%[A-Za-z]');
        if strcmp('Bright',Zone_Type)
            Rec_Bright = load(files{file});
            Rec_Bright = Rec_Bright.Rec_Sigs_B;
        elseif strcmp('Quiet',Zone_Type)
            Rec_Quiet = load(files{file});
            Rec_Quiet = Rec_Quiet.Rec_Sigs_Q;
        else
            error('Error reading the recordings. The type of zone read from the file name is not supported');
        end
        
        if all( [any(Rec_Bright(:)), any(Rec_Quiet(:))] ) % If we have a complete group of speaker signals
            
            ind = strfind(fileName, Zone_Type);
            
            % Read original file
            for file_orig = 1:length(files)
                if ~isempty( strfind(files{file_orig}, [fileName(1:ind-1) 'Original']) )
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
             
            % BEGIN Resize the original speech signal and Align it with the
            % reverberant signals.
            orig(length(orig):size(Rec_Bright,2))=0; % Resize the original signal because the reverberant signal will be longer
            if (length(orig) ~= length(Rec_Bright)) || (length(orig) ~= length(Rec_Quiet))
                error('Size of the original signal does not match the reproduced signal!');
            end
            
            %c_speed = 343;%343m/s speed of sound in air
            %max_delay = speaker_radius*2 / c_speed * Fs;
            max_delay = Fs / 2;
            Original = zeros(Num_Receivers,2,length(orig));
            
            for r = 1:Num_Receivers
                delay = sigalign(Rec_Bright(r,:), orig, [-max_delay 0]) - 1;
                Original(r,1,:) = [orig(-delay:end); zeros(-delay-1,1)];
                
                delay = sigalign( Rec_Quiet(r,:), orig, [-max_delay 0]) - 1;
                Original(r,2,:) = [orig(-delay:end); zeros(-delay-1,1)];
            end
            % END resize and align
            
            
            % Calculate and save results
            
            % Perceptual Evaluation of Speech Quality
            Room_Acoustics.Apply_RIRs.Save_Reverb_PESQ_Result( Original, Rec_Bright, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2) );
            
            % Speech Intelligibility
            Room_Acoustics.Apply_RIRs.Save_Reverb_STOI_Result( Original, Rec_Bright, Rec_Quiet, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2) );
            
            % Signal to Noise Ratio
            Room_Acoustics.Apply_RIRs.Save_Reverb_SNR_Result( Original, Rec_Bright, Rec_Quiet, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2) );
            
            
            
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


