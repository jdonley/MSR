function Save_Reverb_PESQ_Result(Original_, Rec_Sigs_B, Fs, ResultsPath, Output_file_path_ext, FileName)
%SAVE_REVERB_PESQ_RESULT Summary of this function goes here
% 
% Syntax:	SAVE_REVERB_PESQ_RESULT(Original_, Rec_Sigs_B, Rec_Sigs_Q, Fs, ResultsPath, Output_file_path_ext, FileName)
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% 
% 
% See also: 

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015
% Date: 25 August 2015 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_type = 'PESQ';

%% Calculate Results
orig_B = squeeze( Original_(:,1,:) );
for r = 1:size(Rec_Sigs_B,1)
    %PESQ_B(r) = Tools.pesq2( orig_B(r,:), Rec_Sigs_B(r,:), Fs );
    PESQ_B(r) = Tools.pesq3( orig_B(r,:), Rec_Sigs_B(r,:), Fs );
    
    PESQ_MOS_B(r) = pesq2mos( PESQ_B(r) );
end
ConfInt_95 = Tools.confidence_intervals( [PESQ_MOS_B'  PESQ_B'] );

%%
% Work out noise mask
pos = strfind(FileName,'Hz_');
PW_angle = sscanf(FileName(pos:end),'%*[^0-9^-]%d%*s');

pos = strfind(FileName,'Angle');
noise_mask = sscanf(FileName(pos:end),'%*[^0-9^-]%d%*s');

%Determine where to save results
pos = strfind(FileName,'dB');
weight = sscanf(FileName(pos:end),'%*[^0-9]%d%*s');

pos = strfind(FileName,'weight');
mask_type = sscanf(FileName(pos:end),'weight%[^0-9]');
mask_type = ['__' strrep(mask_type,'_','')];

FileName = [FileName(1:pos+5) mask_type];

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

fileID = fopen([ResultsPath Output_file_path_ext results_type '_Results_' num2str(PW_angle) 'pwAngle_' num2str(weight) 'weight' mask_type '.csv'],'a');

fprintf(fileID, ...
        ['%s,' ...
        'Noise_Mask,%f,' ...
        'PESQ_MOS_Bright,%f,' ...
        '95_ConfInt_PESQ_MOS_Bright_L,%f,' ...
        '95_ConfInt_PESQ_MOS_Bright_U,%f,' ...
        'PESQ_Bright,%f,' ...
        '95_ConfInt_PESQ_Bright_L,%f,' ...
        '95_ConfInt_PESQ_Bright_U,%f,' '\r\n'], ...
        FileName, ...
        noise_mask, ...
        mean(PESQ_MOS_B), ...
        ConfInt_95(1,1), ...
        ConfInt_95(1,2), ...
        mean(PESQ_B), ...
        ConfInt_95(2,1), ...
        ConfInt_95(2,2) );
fclose(fileID);

end
