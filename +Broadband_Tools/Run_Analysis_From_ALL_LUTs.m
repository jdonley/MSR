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
%weights = [0 1e0 1e2 1e4];
weights = [10^-0.5 10^0.5 10^1.5 10^2.5];

N_f_low = 16;
N_w_low = 8;
N_increments = 5;
f_lut = repmat(N_f_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1)) ,[N_increments 1]);
w_lut = repmat(N_w_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1))',[1 N_increments]);
f_lut = f_lut(:);
w_lut = w_lut(:);
LUTS = 1:length(f_lut);

for LUT = LUTS
    LUT_tic = tic;
    fprintf('\n====== Building Broadband Multizone Soundfield Zone Signals usings LUT %s ======\n',[num2str(f_lut(LUT)) 'f x ' num2str(w_lut(LUT)) 'w']);
    fprintf('\tCompletion: ');n=0;i=0;
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
            Broadband_Tools.Analyse_Broadband_Signals_from_LUT( ...
                [filePath '\'], ...
                fileName, ...
                fileExt, ...
                weights(w), ...
                f_lut(LUT), ...
                w_lut(LUT));
            
            tElapsed = toc(LUT_tic);
            ratio = (i*length(weights) + w) / (length(weights)*length(files));
            tRem = (1-ratio) / ratio * tElapsed;
            tTot = tElapsed + tRem;
            fprintf(repmat('\b',1,n));
            n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
        end
        i=i+1;
    end
end

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
