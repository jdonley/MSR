function [Freqs,SUPP,ConfInt_Low,ConfInt_Up] = import_Suppression_vs_Freq( SYS )
%import_Suppression_Reverb Import numeric data from a text file as column vectors.
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up]
%   = import_Suppression_Reverb(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up]
%   = import_Suppression_Reverb(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up] = import_Suppression_Reverb('Suppression_Results.mat',1, 180);
%

%%
SYS.analysis_info.f_high = SYS.signal_info.f_high;

SUPP = [];
ConfInt_Low=[];
ConfInt_Up=[];
%%
for a = 1:2
    
    if isfield(SYS.signal_info,'Nfft_optimal') ...
            && ~isempty(SYS.signal_info.Nfft_optimal)
        SYS.signal_info.Nfft = SYS.signal_info.Nfft_optimal;
    else
        % Assume 12ms is optimal following from:
        % J. Donley, C. Ritz & W. B. Kleijn, "Active Speech Control using
        % Wave-Domain Processing with a Linear Wall of Dipole Secondary
        % Sources," in International Conference on Acoustics, Speech and
        % Signal Processing (ICASSP). IEEE, 2017, pp. 456-460.
        SYS.signal_info.Nfft = 12 * 1e-3 * SYS.signal_info.Fs;
    end
    
    if a==1
        if isfield(SYS.signal_info,'system_time_delay') ...
                && ~isempty(SYS.signal_info.system_time_delay)
            SYS.signal_info.time_delay = SYS.signal_info.system_time_delay;
        else
            SYS.signal_info.time_delay = []; %If empty the time delay will based on the frame length
        end
    elseif a==2
        SYS.signal_info.time_delay = [0]; %If empty the time delay will based on the frame length
    end
    SYS.signal_info.zeropadtime = SYS.signal_info.Nfft / SYS.signal_info.Fs; % miliseconds of zero padding (i.e. maximum time shift before circular convolution)
    SYS.system_info.LUT_frequencies = (SYS.signal_info.Nfft + SYS.signal_info.zeropadtime * SYS.signal_info.Fs)/2;
    SYS.system_info.LUT_resolution = [num2str(SYS.system_info.LUT_frequencies) 'f' ...
        SYS.system_info.sc ...
        num2str(SYS.system_info.LUT_weights) 'w'];
    
    %%
    SYS.Main_Setup = SYS.Main_Setup(1);
%         Falias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c;
%     Falias = 2000;
    
    %%
    SYS.signal_info.method = SYS.signal_info.methods_list{end};
    load([Results.getResultsPath(SYS) 'Suppression_Results.mat'])
    A = [S{:}]; sk = size(S{1},2);
    B = [A{2:sk:end}];
    C = [A{3:sk:end}];
    C2=[];
    for i=1:20
        tmp= Tools.confidence_intervals( [A{sk*(i-1)+1}].' , 95);
        C2(:,i)= tmp(:,2);
    end
    D = [A{sk:sk:end}];
    
    SUPP(:,a) = mean(B,2).';
    ConfInt_Low(:,a) = mean(C2,2).';
    
    
end

f_ = S{1}{4};
f_(1)=[];SUPP(1,:)=[];ConfInt_Low(1,:)=[];
fi=(f_>=SYS.analysis_info.f_low & f_<=SYS.analysis_info.f_high);
Freqs=f_(fi).' / 1e3; %KiloHertz
SUPP(~fi,:)=[];
ConfInt_Low(~fi,:)=[];

ConfInt_Up = ConfInt_Low;



