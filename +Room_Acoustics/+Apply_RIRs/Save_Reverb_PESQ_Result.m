function Save_Reverb_PESQ_Result(Original_, Rec_Sigs_B, Fs, noise_mask, ResultsPath, Output_file_path_ext, FileName, pesqNumber)
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
if nargin < 7
    pesqNumber = 0;
end
results_type = 'PESQ';

%% Calculate Results
orig_B = permute( Original_(:,1,:) ,[1 3 2]);
parfor r = 1:size(Rec_Sigs_B,1)
    PESQ_MOS_B(r) = Tools.pesq_mex_vec( orig_B(r,:), Rec_Sigs_B(r,:), Fs );
%     PESQ_MOS_B(r) = Tools.pesq_mex_fast_vec( orig_B(r,:), Rec_Sigs_B(r,:), Fs );
end
ConfInt_95 = Tools.confidence_intervals( [PESQ_MOS_B'] );

%% Determine where to save results

%[ PW_angle, noise_mask, weight, mask_type ] = Results.getInfoFromFilename(FileName);

%pos = strfind(FileName,'weight');
%FileName = [FileName(1:pos+5) mask_type];

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

fileID = fopen([ResultsPath Output_file_path_ext results_type '_Results.csv'],'a');

fprintf(fileID, ...
        ['%s,' ...
        'Noise_Mask,%f,' ...
        'PESQ_MOS_Bright,%f,' ...
        '95_ConfInt_PESQ_MOS_Bright_L,%f,' ...
        '95_ConfInt_PESQ_MOS_Bright_U,%f,' '\r\n'], ...
        FileName, ...
        noise_mask, ...
        mean(PESQ_MOS_B), ...
        ConfInt_95(1,1), ...
        ConfInt_95(1,2) );
fclose(fileID);

end
