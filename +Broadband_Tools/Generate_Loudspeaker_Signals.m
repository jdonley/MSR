%% Initialise
clc;
clear;
close all;
tic;
% Kill dropbox
%Tools.Dropbox('kill');
%if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
%matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))


%% Analyse broadband signals script
LUT_res = '512f_256w';

loudspeakers = 295;
planewave_angle = 0;

Input_file_path = '+Miscellaneous\+Speech_Files\';
files = Tools.getAllFiles(Input_file_path);
%weights = [0 1e0 1e2 1e4];
weights = [1e4];
noise_masks = [-40 -35 -30 -25 -20 -15 -10 -5 -0];

%mask_type = 'ZoneWeightMask';
mask_type = 'FlatMask';
%mask_type = 'NoMask';




%%
fprintf('\n====== Building Broadband Multizone-Soundfield Loudspeaker-Signals ======\n');
fprintf(['Privacy Weighting: ' mask_type '\n']);n=0;
fprintf('\tCompletion: ');
%parfor_progress( length(files) * length(weights) * length(noise_masks));
for file = 1:length(files)
    try
        audioread(files{file});
    catch err
        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
            continue; % Skip unsupported files
        end
    end
    
    
    for w = 1:length(weights)
        [filePath, fileName, fileExt] = fileparts(files{file});
        
        if strcmp(mask_type, 'NoMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT( ...
                    [filePath '\'], ...
                    fileName, ...
                    fileExt, ...
                    LUT_res, ...
                    weights(w), ...
                    loudspeakers, ...
                    planewave_angle);
        end
        for m = 1:length(noise_masks)
            
            if strcmp(mask_type, 'FlatMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_FlatMask( ...
                    [filePath '\'], ...
                    fileName, ...
                    fileExt, ...
                    LUT_res, ...
                    noise_masks(m), ...
                    weights(w), ...
                    loudspeakers, ...
                    planewave_angle); %Planewave Angle in Bright Zone
                
            elseif strcmp(mask_type, 'ZoneWeightMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask( ...
                    [filePath '\'], ...
                    fileName, ...
                    fileExt, ...
                    LUT_res, ...
                    noise_masks(m), ...
                    weights(w), ...
                    loudspeakers, ...
                    planewave_angle); %Planewave Angle in Bright Zone
            end
            
            %parfor_progress;
        end
        
    end
    
    
end
%parfor_progress(0);

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
%Tools.Dropbox('start');