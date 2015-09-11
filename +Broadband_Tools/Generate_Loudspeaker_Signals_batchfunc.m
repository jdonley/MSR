function Generate_Loudspeaker_Signals_batchfunc(Zone_Weights, Noise_Mask_Levels, Masker_Type, Planewave_Angle)

%% Initialise
% clc;
% clear;
% close all;
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

Input_file_path = '+Miscellaneous\+Speech_Files\';
%Input_file_path = '+Miscellaneous\+TestAudio_Files\';
files = Tools.getAllFiles(Input_file_path);


%%
fprintf('\n====== Building Broadband Multizone-Soundfield Loudspeaker-Signals ======\n');
fprintf(['         Zone Weights: ' strrep(sprintf(strrep(repmat('%d',1,length(Zone_Weights)),'d%','d %'),Zone_Weights),' ',' & ') '\n']);
fprintf(['    Noise Mask Levels: ' [strrep(sprintf(strrep(repmat('%d',1,length(Noise_Mask_Levels)),'d%','d %'),Noise_Mask_Levels),' ','dB & ') 'dB'] '\n']);
fprintf(['          Masker Type: ' Masker_Type '\n']);
fprintf(['      Planewave Angle: ' num2str(Planewave_Angle) '\n\n']);n=0;
fprintf('\tCompletion: ');

F=length(files);
for file = 1:F
    try
        audioread(files{file});
    catch err
        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
            continue; % Skip unsupported files
        end
    end
    
    W = length(Zone_Weights);
    for w = 1:W
        [filePath, fileName, fileExt] = fileparts(files{file});
        
        if strcmp(Masker_Type, 'NoMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT( ...
                    [filePath '\'], ...
                    fileName, ...
                    fileExt, ...
                    LUT_res, ...
                    Zone_Weights(w), ...
                    loudspeakers, ...
                    Planewave_Angle);
        end
        M = length(Noise_Mask_Levels);
        for m = 1:M
            
            if strcmp(Masker_Type, 'FlatMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_FlatMask( ...
                    [filePath '\'], ...
                    fileName, ...
                    fileExt, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    loudspeakers, ...
                    Planewave_Angle); %Planewave Angle in Bright Zone
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask( ...
                    [filePath '\'], ...
                    fileName, ...
                    fileExt, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    loudspeakers, ...
                    Planewave_Angle); %Planewave Angle in Bright Zone
            end
            
            n = Tools.showTimeToCompletion( ((file-1)*W*M + (w-1)*M + m) / (F*W*M), n);
        end
        
    end
    
    
end


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
%Tools.Dropbox('start');