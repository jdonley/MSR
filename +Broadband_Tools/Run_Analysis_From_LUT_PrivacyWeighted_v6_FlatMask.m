%% Initialise
clc;
clear;
close all;
tic;
%if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
%matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))


%% Analyse broadband signals script
Input_file_path = '+Miscellaneous\+Speech_Files\';
files = Tools.getAllFiles(Input_file_path);

fprintf('\n====== Building Broadband Multizone-Soundfield Zone-Signals ======\n');
fprintf(['Privacy Weighting: ' 'v6 with Flat Mask' '\n']);n=0;
fprintf('\tCompletion: ');n=0;
Noise_Masks = 0:5:45;

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
        Broadband_Tools.Broadband_PrivacyWeighted_v6_FlatMask( ...
            [filePath '\'], ...
            fileName, ...
            fileExt, ...
            '512f_256w', ...
            Noise_Masks(m), ... 
            295); %295 loudspeaker setup
            %65); %65 loudspeaker setup
    end
    tElapsed = toc;
    ratio = m/length(Noise_Masks);
    tRem = (1-ratio) / ratio * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
end

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
