function Save_Reverb_Suppression_Result(Talker_, Rec_Sigs_B, Fs, Nfft, f_band, ResultsPath, Output_file_path_ext, FileName)
%SAVE_REVERB_SUPPRESSION_RESULT Summary of this function goes here
% 
% Syntax:	SAVE_REVERB_SUPPRESSION_RESULT(Talker_, Rec_Sigs_B, Rec_Sigs_Q, Fs, Nfft, ResultsPath, Output_file_path_ext, FileName)
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
% Copyright: Jacob Donley 2016
% Date: 8 August 2016 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_type = 'Suppression';

%% Calculate Results
% T = fft(Talker_.',Nfft);
% T = abs(T(1:end/2+1));
Tpow = pwelch(Talker_.',hamming(Nfft),Nfft/2,Nfft,Fs,'power');
T = sqrt(Tpow);
% C = fft(Rec_Sigs_B.',Nfft);
% C = abs(C(1:end/2+1));
Cpow = pwelch(Rec_Sigs_B.',hamming(Nfft),Nfft/2,Nfft,Fs,'power');
C = sqrt(Cpow);

SUPP_B = mag2db(C) - mag2db(T);
frqs = linspace(0,1,Nfft/2+1)*Fs/2;

%%
Supp_Mean = mean(SUPP_B,2);
Supp_Std = sqrt( var(SUPP_B,0,2) );

%%
f_band_indices = frqs>=f_band(1) & frqs <= f_band(2);
Supp_BandMean = mean(Supp_Mean(f_band_indices));
Supp_BandStd = mean(Supp_Std(f_band_indices));

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

fn = [ResultsPath Output_file_path_ext results_type '_Results'];
if ~exist([fn '.mat'],'file'); S={};save(fn,'S'); end

D = load([ResultsPath Output_file_path_ext results_type '_Results']);S=D.S;
S{end+1} = {...
    SUPP_B,...
    Supp_Mean,...
    Supp_Std,...
    frqs,...
    FileName,...
    {Supp_BandMean, Supp_BandStd, f_band}   };
save([ResultsPath Output_file_path_ext results_type '_Results'],'S');

% fileID = fopen([ResultsPath Output_file_path_ext results_type '_Results.csv'],'a');

% fprintf(fileID, ...
%         ['%s,' ...
%         'Noise_Mask,%f,' ...
%         'PESQ_MOS_Bright,%f,' ...
%         '95_ConfInt_PESQ_MOS_Bright_L,%f,' ...
%         '95_ConfInt_PESQ_MOS_Bright_U,%f,' '\r\n'], ...
%         FileName, ...
%         noise_mask, ...
%         mean(SUPP_B), ...
%         ConfInt_95(1,1), ...
%         ConfInt_95(1,2) );
% fclose(fileID);

end
