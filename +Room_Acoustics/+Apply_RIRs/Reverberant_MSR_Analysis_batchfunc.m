function Reverberant_MSR_Analysis_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, pw_angle, mask_level, pesqNumber)
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

%loudspeakers = 295; %Number of loudspeakers
%loudspeakers = 32; %Number of loudspeakers
loudspeakers = 24; %Number of loudspeakers
%loudspeakers = 16; %Number of loudspeakers

loudspeaker_layout = {'numberof_loudspeakers',        loudspeakers, ...
                      'loudspeaker_radius',           1.5, ...
                      'angleof_loudspeakerarc'        180, ...
                      'loudspeaker_model',            'Genelec 8010A', ...
                      'angleof_loudspeakerarrcentre', 180, ...
                      'loudspeaker_spacing',          0.01    };            

zone_layout = {'brightzone_source_angle',     pw_angle};
             
setup = Speaker_Setup.createSetup({...
            zone_layout{:}, ...
            loudspeaker_layout{:}});
        
        
Num_Receivers = 32; %Number of recording points (microphones)

Reproduction_Centre = Room_Size ./ 2;

if nargin < 5
   mask_level = [];
end
if nargin < 6
    pesqNumber = 0;
end

Output_file_path_ext = ['+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc_LUT_' LUT_resolution '\'];

%% Load RIR Database file
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

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


%% Delete previous results file
results_type{1} = 'PESQ';
results_type{2} = 'SI';
results_type{3} = 'SNR';

[~, FileName, ~] = fileparts(files{1});
Num_ = strfind(FileName,'__')+2;
Num_p = Num_(3);
z_t = sscanf([' ' FileName(Num_p:end)],'%*[^_]_%[A-Za-z]');
ind = strfind(FileName, z_t);
[ pw_a, ~, wei, m_t ] = Results.getInfoFromFilename(FileName(1:ind-2));

for i=1:length(results_type)
    csvFilename = ([ResultsPath Output_file_path_ext results_type{i} '_Results_' num2str(pw_a) 'pwAngle_' num2str(wei) 'weight' m_t '.csv']);
    delete(csvFilename);
end

%% Start Evaluation Loop

fprintf('\n====== Analysing Simulated Reverberant Signals ======\n');
fprintf(['            Room Size: ' [strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','m x ') 'm'] '\n']);
fprintf(['Wall Absorption Coeff: ' num2str(Wall_Absorption_Coeff) '\n']);
fprintf(['      Planewave Angle: ' num2str(pw_angle) '\n']);
fprintf(['    Privacy Weighting: ' mask_type '\n\n']);n=0;
fprintf('\tCompletion: ');


fileName_prev = '';
Rec_Bright = [];
Rec_Quiet = [];
for file = 1:length(files)
   
    [~, fileName, ~] = fileparts(files{file});
    
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
             
            % BEGIN Resize the original speech signal and
            % Align it with the reverberant signals.
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
            
            % BEGIN Calculate and save results
            
            % Perceptual Evaluation of Speech Quality
            Room_Acoustics.Apply_RIRs.Save_Reverb_PESQ_Result( Original, Rec_Bright, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2), pesqNumber );
            
            % Speech Intelligibility
                % STOI
            Room_Acoustics.Apply_RIRs.Save_Reverb_STOI_Result( Original, Rec_Bright, Rec_Quiet, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2) );
                % STI
            %Room_Acoustics.Apply_RIRs.Save_Reverb_STI_Result( Original, Rec_Bright, Rec_Quiet, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2) );
            
            % Signal to Noise Ratio
            Room_Acoustics.Apply_RIRs.Save_Reverb_SNR_Result( Original, Rec_Bright, Rec_Quiet, Fs, ResultsPath, Output_file_path_ext, fileName(1:ind-2) );
            
            % END calc and save results
            
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

