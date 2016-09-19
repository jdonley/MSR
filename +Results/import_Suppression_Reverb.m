function [Noise_Mask_Level,PESQ_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] = import_Suppression_Reverb(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [DOMAIN_VEC,SUPP_BRIGHT,CONFINT_BRIGHT_LOW,CONFINT_BRIGHT_UP]
%   = import_Suppression_Reverb(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [NOISE_MASK_LEVEL,PESQ_BRIGHT,CONFINT_BRIGHT_LOW,CONFINT_BRIGHT_UP]
%   = import_PESQ_Reverb(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [Noise_Mask_Level,PESQ_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] = import_PESQ_Reverb('PESQ_Results_10000weight__withFlatMask.csv',1, 180);
%
%    See also TEXTSCAN.


%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column3: double (%f)
%	column5: double (%f)
%   column7: double (%f)
%	column9: double (%f)
%   column11: double (%f)
%	column13: double (%f)
%   column15: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%*s%f%*s%f%*s%f%[^\n\r]';

%% Open the text file.
[fileID, errMSG] = fopen(filename,'r');
if fileID ~= -1
    %% Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    %% Close the text file.
    fclose(fileID);
    
    %% Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.
    
    %% Allocate imported array to column variable names
    Noise_Mask_Level = dataArray{:, 1};
    PESQ_Bright = dataArray{:, 2};
    ConfInt_Bright_Low = dataArray{:, 3};
    ConfInt_Bright_Up = dataArray{:, 4};
    
else
    Noise_Mask_Level = fileID;
    PESQ_Bright = errMSG;
    ConfInt_Bright_Low = [];
    ConfInt_Bright_Up = [];
end
