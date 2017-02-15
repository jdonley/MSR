function Save_Reverb_SPLs_Result(Rec_Sigs_B, Rec_Sigs_Q, Fs, Nfft, f_band, ResultsPath, Output_file_path_ext, FileName, SYS)
%SAVE_REVERB_SPLS_RESULT Summary of this function goes here
% 
% Syntax:	SAVE_REVERB_SPLS_RESULT(Rec_Sigs_B, Rec_Sigs_Q, Fs, Nfft, ResultsPath, Output_file_path_ext, FileName, SYS)
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
results_type = 'SPL';

%% Deconvolve responses
invFilt = load(cell2mat(Tools.getAllFiles(SYS.signal_info.InverseFilter_filepath)));invFilt=invFilt.invY;
irQ = Tools.extractIR(Rec_Sigs_Q,invFilt);
irB = Tools.extractIR(Rec_Sigs_B,invFilt);

%% Calculate Results
% T = fft(Talker_.',Nfft);
% T = abs(T(1:end/2+1));
Qpow = pwelch(irQ,hamming(Nfft),Nfft/2,Nfft,Fs,'power');
Q = sqrt(Qpow);
% C = fft(Rec_Sigs_B.',Nfft);
% C = abs(C(1:end/2+1));
Bpow = pwelch(irB,hamming(Nfft),Nfft/2,Nfft,Fs,'power');
B = sqrt(Bpow);

SPL_B = mag2db(B);
SPL_Q = mag2db(Q);
frqs = linspace(0,1,Nfft/2+1)*Fs/2;

%%
SPL_B_Mean = mean(SPL_B,2);
SPL_B_Std = sqrt( var(SPL_B,0,2) );
SPL_Q_Mean = mean(SPL_Q,2);
SPL_Q_Std = sqrt( var(SPL_Q,0,2) );

%%
f_band_indices = frqs>=f_band(1) & frqs <= f_band(2);
SPLB_BandMean = mean(SPL_B_Mean(f_band_indices));
SPLB_BandStd = mean(SPL_B_Std(f_band_indices));
SPLQ_BandMean = mean(SPL_Q_Mean(f_band_indices));
SPLQ_BandStd = mean(SPL_Q_Std(f_band_indices));

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

fn = [ResultsPath Output_file_path_ext results_type '_Results'];
if ~exist([fn '.mat'],'file'); S={};save(fn,'S'); end

D = load([ResultsPath Output_file_path_ext results_type '_Results']);S=D.S;
S{end+1} = {...
    SPL_B,...
    SPL_B_Mean,...
    SPL_B_Std,...
    SPL_Q,...
    SPL_Q_Mean,...
    SPL_Q_Std,...
    frqs,...
    FileName,...
    {SPLB_BandMean, SPLB_BandStd, SPLQ_BandMean, SPLQ_BandStd, f_band}   };
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
