function Save_Reverb_STOI_Result(Original_, Rec_Sigs_B, Rec_Sigs_Q, Fs, ResultsPath, Output_file_path_ext, FileName)
%SAVE_REVERB_STOI_RESULT Summary of this function goes here
% 
% Syntax:	SAVE_REVERB_STOI_RESULT(Original_, Rec_Sigs_B, Rec_Sigs_Q, Fs, ResultsPath, Output_file_path_ext, FileName)
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% 
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015
% Date: 05 August 2015 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_type = 'STI_IR';

%% Calculate Results
orig_B = squeeze( Original_(:,1,:) );
orig_Q = squeeze( Original_(:,2,:) );

invsweepfft = load('+Tools\invsweepfft.mat', 'invsweepfft'); invsweepfft = invsweepfft.invsweepfft;
parfor r = 1:size(Rec_Sigs_Q,1)
    [IR, IRinv] = Tools.extractIR(Rec_Sigs_B(r,:), invsweepfft);
    ImpResp_Bright(r,:) = [IRinv; IR];
    [IR, IRinv] = Tools.extractIR(Rec_Sigs_Q(r,:), invsweepfft);
    ImpResp_Quiet(r,:) = [IRinv; IR];
end

%ImpResp_Bright = mean(ImpResp_Bright,1);
%ImpResp_Quiet = mean(ImpResp_Quiet,1);

maxVal = max(abs( [ImpResp_Bright(:); ImpResp_Quiet(:)] ));
ImpResp_Bright = ImpResp_Bright / maxVal;
ImpResp_Quiet = ImpResp_Quiet / maxVal;

%% Determine where to save results
% Work out Planewave Angle
pos = strfind(FileName,'Hz_');
PW_angle = sscanf(FileName(pos:end),'%*[^0-9^-]%d%*s');

if ~isempty(strfind(FileName,'with'))
    % Work out noise mask
    pos = strfind(FileName,'Angle');
    noise_mask = sscanf(FileName(pos:end),'%*[^0-9^-]%d%*s');
    
    % Work out weight
    pos = strfind(FileName,'dB');
    weight = sscanf(FileName(pos:end),'%*[^0-9]%d%*s');
    
    % Work out mask type
    pos = strfind(FileName,'weight');
    mask_type = sscanf(FileName(pos:end),'weight%[^0-9]');
    mask_type = ['__' strrep(mask_type,'_','')];
else
    % Work out weight
    pos = strfind(FileName,'Angle');
    weight = sscanf(FileName(pos:end),'%*[^0-9]%d%*s');
    
    % Assign mask type and noise mask
    mask_type = 'NoMask';
    noise_mask = -Inf;
end

%% Aave Speech Transmission Index Impulse Repsonse to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

for r = 1:1:size(Rec_Sigs_Q,1)
    filepath_B = [ResultsPath Output_file_path_ext results_type '_Results_' num2str(PW_angle) 'pwAngle_' num2str(noise_mask) 'dB_' num2str(weight) 'weight' mask_type '__Bright_mic' num2str(r) '.WAV'];
    filepath_Q = [ResultsPath Output_file_path_ext results_type '_Results_' num2str(PW_angle) 'pwAngle_' num2str(noise_mask) 'dB_' num2str(weight) 'weight' mask_type '__Quiet_mic' num2str(r) '.WAV'];
    
    audiowrite(filepath_B, ImpResp_Bright(r,:), Fs);
    audiowrite(filepath_Q, ImpResp_Quiet(r,:), Fs);
end

end
