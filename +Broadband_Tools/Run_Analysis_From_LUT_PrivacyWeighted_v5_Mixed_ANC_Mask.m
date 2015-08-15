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
Input_file_ext  = '.wav';
files = Tools.getAllFiles(Input_file_path);

fprintf('\n====== Building Broadband Multizone-Soundfield Zone-Signals ======\n');
fprintf(['Privacy Weighting: ' 'v5_Mixed_ANC_Mask' '\n']);
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
        Broadband_Tools.Broadband_PrivacyWeighted_v5_Mixed_ANC_Mask( ...
            [filePath '\'], ...
            fileName, ...
            fileExt, ...
            '512f_256w', ... % LUT resolution
            Noise_Masks(m), ... 
            65); %65 loudspeaker setup
        
        
        tElapsed = toc;
        ratio = ((m-1)*length(files)+file)/(length(Noise_Masks)*length(files));
        tRem = (1-ratio) / ratio * tElapsed;
        tTot = tElapsed + tRem;
        fprintf(repmat('\b',1,n));
        n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
    end
end

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
Tools.Dropbox('start');