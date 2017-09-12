function Save_Reverb_RIR_Result(Talker, Rec_Sigs, ResultsPath, Output_file_path_ext, FileName, SYS)
%SAVE_REVERB_RIR_RESULT Summary of this function goes here
% 
% Syntax:	SAVE_REVERB_RIR_RESULT(Rec_Sigs_B, Rec_Sigs_Q, Fs, Nfft, ResultsPath, Output_file_path_ext, FileName, SYS)
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
% Copyright: Jacob Donley 2017
% Date: 10 September 2017
% Revision: 0.1 (10 September 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_type = 'RIR';

%% Filter to analysis frequency range
fband = [SYS.analysis_info.f_low SYS.analysis_info.f_high];
[b,a] = cheby1(5,1,fband/SYS.signal_info.Fs);
T = filter(b,a,Talker);
R = filter(b,a,Rec_Sigs);

%% Deconvolve responses
invFilt = load(cell2mat(Tools.getAllFiles(SYS.signal_info.InverseFilter_filepath)));invFilt=invFilt.invY;
irT = Tools.extractIR(T,invFilt);
irR = Tools.extractIR(R,invFilt);

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist([ResultsPath Output_file_path_ext],'dir'); mkdir([ResultsPath Output_file_path_ext]); end

fn = [ResultsPath Output_file_path_ext results_type '_Results'];
if ~exist([fn '.mat'],'file'); S={};save(fn,'S'); end

D = load([ResultsPath Output_file_path_ext results_type '_Results']);S=D.S;
S{end+1} = {...
    irT,...
    irR,...
    FileName };
save(fn,'S');


end
