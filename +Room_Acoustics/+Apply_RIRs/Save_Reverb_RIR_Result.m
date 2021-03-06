 function Save_Reverb_RIR_Result(Talker, TalkerAnechoic, Rec_Sigs, ResultsPath, Output_file_path_ext, FileName, SYS)
%SAVE_REVERB_RIR_RESULT Summary of this function goes here
% 
% Syntax:	Save_Reverb_RIR_Result(Talker, TalkerAnechoic, Rec_Sigs, ResultsPath, Output_file_path_ext, FileName, SYS)
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

%% Get the cancellation signal only
Cancel_Sig = Rec_Sigs - Talker;

%% Get Reverberant part only
Talker_Reverb = Talker - TalkerAnechoic;
Rec_Sigs_Reverb = Rec_Sigs - TalkerAnechoic;

%% Filter to analysis frequency range
 fband = [SYS.analysis_info.f_low SYS.analysis_info.f_high];
[b,a] = cheby1(5,1,fband/(SYS.signal_info.Fs/2));

Cancel_Sig = filter(b,a,Cancel_Sig);

Talker = filter(b,a,Talker);
Rec_Sigs = filter(b,a,Rec_Sigs);

Talker_Reverb = filter(b,a,Talker_Reverb);
Rec_Sigs_Reverb = filter(b,a,Rec_Sigs_Reverb);


%% Deconvolve responses
invFilt = load(cell2mat(Tools.getAllFiles(SYS.signal_info.InverseFilter_filepath)));invFilt=invFilt.invY;

irC = Tools.extractIR(Cancel_Sig,invFilt);

irT = Tools.extractIR(Talker,invFilt);
irR = Tools.extractIR(Rec_Sigs,invFilt);

irTrvrb = Tools.extractIR(Talker_Reverb,invFilt);
irRrvrb = Tools.extractIR(Rec_Sigs_Reverb,invFilt);


%%%%
Fs = SYS.signal_info.Fs;
[p1,ff]=pwelch( (imag(hilbert(irTrvrb))),hamming(1024,'p'),512,1024,Fs,'power');
p2 = pwelch( (imag(hilbert(irTrvrb)) + irC*2.899),hamming(1024,'p'),512,1024,Fs,'power');
figure(1);
plot(ff/1e3,pow2db([p1,p2]));
set(gca,'xscale','log');grid on;
xlim(fband/1e3);
figure(2);
plot(linspace(0,Fs,numel(Talker_Reverb))/1e3,...
    (mod((angle(fft(Talker_Reverb))) - (angle(fft(Cancel_Sig))),2*pi))/pi*180-180);
set(gca,'xscale','log');grid on;
xlim(fband/1e3);
hold on
plot([0.1 8],[90 90]); hold off;
%%%%

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
