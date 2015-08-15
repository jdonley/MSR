clc;
clear;
close all;

%% Analyse broadband signals script
Input_file_path = '+Miscellaneous\+Speech_Files\';
files = Tools.getAllFiles(Input_file_path);
%weights = [0 1e0 1e2 1e4];
weights = [10^-0.5 10^0.5 10^1.5 10^2.5];

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
         Broadband_Tools.Analyse_Broadband_Signals( ...
            [filePath '\'], ...
            fileName, ...
            fileExt, ...
            weights(w));
    end
end

fprintf('\n\n FINISHED \n\n');