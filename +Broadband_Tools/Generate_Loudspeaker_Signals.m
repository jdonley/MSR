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
Input_file_path = '+Miscellaneous\+Speech_Files\';
files = Tools.getAllFiles(Input_file_path);
weights = [0 1e0 1e2 1e4];
weights = [1e4];
noise_masks = [-40 -35 -30 -25 -20 -15 -10 -5 -0];

fprintf('\n====== Building Broadband Multizone-Soundfield Loudspeaker-Signals ======\n');
fprintf(['Privacy Weighting: ' 'None (Uniform Weighting)' '\n']);n=0;
fprintf('\tCompletion: ');
parfor_progress( length(files) * length(weights) * length(noise_masks));
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
        
%         Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT( ...
%             [filePath '\'], ...
%             fileName, ...
%             fileExt, ...
%             '512f_256w', ...
%             weights(w), ...
%             295);
        
        for m = 1:length(noise_masks)
            Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_FlatMask( ...
                [filePath '\'], ...
                fileName, ...
                fileExt, ...
                '512f_256w', ...
                noise_masks(m), ...
                weights(w), ...
                295, ...
                15); %Planewave Angle in Bright Zone

%         for m = 1:length(noise_masks)
%             Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask( ...
%                 [filePath '\'], ...
%                 fileName, ...
%                 fileExt, ...
%                 '512f_256w', ...
%                 noise_masks(m), ...
%                 weights(w), ...
%                 295, ...
%                 15); %Planewave Angle in Bright Zone
            
            parfor_progress;
        end
        
    end
    
    
    
    
    
    
    
end
parfor_progress(0);

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
%Tools.Dropbox('start');