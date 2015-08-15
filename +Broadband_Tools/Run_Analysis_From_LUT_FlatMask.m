%% Initialise
clc;
clear;
close all;
tic;
% Kill dropbox
Tools.Dropbox('kill');
%if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
%matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))


%% Analyse broadband signals script
Input_file_path = '+Miscellaneous\+Speech_Files\';
files = Tools.getAllFiles(Input_file_path);
weights = [0 1e0 1e2 1e4];

fprintf('\n====== Building Broadband Multizone-Soundfield Zone-Signals ======\n');
fprintf(['Privacy Weighting: ' 'None (Uniform Weighting)' '\n']);n=0;
fprintf('\tCompletion: ');

Noise_Masks = 0:5:45;
%weight = 0;
weight = 1e4;

parfor_progress( length(files) * length(Noise_Masks) );

for m = 1:length(Noise_Masks)
    for file = 1:length(files)
        try
            audioread(files{file});
        catch err
            if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
                continue; % Skip unsupported files
            end
        end
        
        [filePath, fileName, fileExt] = fileparts(files{file});
        Broadband_Tools.Analyse_Broadband_Signals_from_LUT_FlatMask( ...
            [filePath '\'], ...
            fileName, ...
            fileExt, ...
            '512f_256w', ...
            Noise_Masks(m), ...
            weight, ... 
            65); %65 loudspeaker setup
        parfor_progress;
    end
end
parfor_progress(0);

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
Tools.Dropbox('start');