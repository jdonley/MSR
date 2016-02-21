function Save_Reverb_SNR_Result(Original_, Rec_Sigs_B, Rec_Sigs_Q, Fs, noise_mask, ResultsPath, Output_file_path_ext, FileName)
%SAVE_REVERB_SNR_RESULT Summary of this function goes here
% 
% Syntax:	SAVE_REVERB_SNR_RESULT(Original_, Rec_Sigs_B, Rec_Sigs_Q, Fs, noise_mask, ResultsPath, Output_file_path_ext, FileName)
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
% Date: 24 August 2015 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_type = 'SNR';

%% Calculate Results
orig_B = squeeze( Original_(:,1,:) );
orig_Q = squeeze( Original_(:,2,:) );
parfor r = 1:size(Rec_Sigs_Q,1)
    [SNRSeg_B(r), SNR_B(r)] =  snrseg( Rec_Sigs_B(r,:), orig_B(r,:), Fs );
    [SNRSeg_Q(r), SNR_Q(r)] =  snrseg( Rec_Sigs_Q(r,:), orig_Q(r,:), Fs );
end
ConfInt_95 = Tools.confidence_intervals( [SNRSeg_B' SNRSeg_Q'  SNR_B' SNR_Q'] );

%% Determine where to save results
   
% Work out Planewave Angle
% pos = strfind(FileName,'Hz_');
% PW_angle = sscanf(FileName(pos:end),'%*[^0-9^-]%d%*s');
% 
% if ~isempty(strfind(FileName,'with'))
%     % Work out noise mask
%     pos = strfind(FileName,'Angle');
%     noise_mask = sscanf(FileName(pos:end),'%*[^0-9^-]%d%*s');
%     
%     % Work out weight
%     pos = strfind(FileName,'dB');
%     weight = sscanf(FileName(pos:end),'%*[^0-9]%d%*s');
%     
%     % Work out mask type
%     pos = strfind(FileName,'weight');
%     mask_type = sscanf(FileName(pos:end),'weight%[^0-9]');
%     mask_type = ['__' strrep(mask_type,'_','')];
% else
%     % Work out weight
%     pos = strfind(FileName,'Angle');
%     weight = sscanf(FileName(pos:end),'%*[^0-9]%d%*s');
%     
%     % Assign mask type and noise mask
%     mask_type = 'NoMask';
%     noise_mask = -Inf;
% end
% 
% pos = strfind(FileName,'weight');
% FileName = [FileName(1:pos+5) mask_type];

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

fileID = fopen([ResultsPath Output_file_path_ext results_type '_Results.csv'],'a');

fprintf(fileID, ...
        ['%s,' ...
        'Noise_Mask,%f,' ...
        'SNRSeg_dB_Bright,%f,' ...
        'SNRSeg_dB_Quiet,%f,' ...
        '95_ConfInt_SNRSeg_Bright_L,%f,' ...
        '95_ConfInt_SNRSeg_Bright_U,%f,' ...
        '95_ConfInt_SNRSeg_Quiet_L,%f,' ...
        '95_ConfInt_SNRSeg_Quiet_U,%f,' ...
        'SNR_dB_Bright,%f,' ...
        'SNR_dB_Quiet,%f,' ...
        '95_ConfInt_SNR_Bright_L,%f,' ...
        '95_ConfInt_SNR_Bright_U,%f,' ...
        '95_ConfInt_SNR_Quiet_L,%f,' ...
        '95_ConfInt_SNR_Quiet_U,%f,'  '\r\n'], ...
        FileName, ...
        noise_mask, ...
        mean(SNRSeg_B), ...
        mean(SNRSeg_Q), ...
        ConfInt_95(1,1), ...
        ConfInt_95(1,2), ...
        ConfInt_95(2,1), ...
        ConfInt_95(2,2), ...
        mean(SNR_B), ...
        mean(SNR_Q), ...
        ConfInt_95(3,1), ...
        ConfInt_95(3,2), ...
        ConfInt_95(4,1), ...
        ConfInt_95(4,2) );
fclose(fileID);

end
