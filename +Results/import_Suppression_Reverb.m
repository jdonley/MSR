function [BlockLen,SUPP,ConfInt_Low,ConfInt_Up] = import_Suppression_Reverb( SYS )
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


N_vec = (4:4:32).'*1e-3*SYS.signal_info.Fs;
td_vec = {0,[]};

SUPP = zeros(numel(N_vec),numel(td_vec));
ConfInt_Low = SUPP;

for td_ = 1:numel(td_vec)
    for N_ = 1:numel(N_vec)
        SYS.signal_info.Nfft = N_vec(N_);
        SYS.signal_info.time_delay = td_vec{td_};
        SYS.signal_info.zeropadtime = SYS.signal_info.Nfft / SYS.signal_info.Fs;
        SYS.system_info.LUT_frequencies = (SYS.signal_info.Nfft + SYS.signal_info.zeropadtime * SYS.signal_info.Fs)/2;
        SYS.system_info.LUT_resolution = [num2str(SYS.system_info.LUT_frequencies) 'f' ...
                              SYS.system_info.sc ...
                              num2str(SYS.system_info.LUT_weights) 'w'];
        SYS.Main_Setup = SYS.Main_Setup(1);
        SYS.signal_info.method = SYS.signal_info.methods_list{end};
        load([Results.getResultsPath(SYS) 'Suppression_Results.mat']);
        A = [S{:}]; sk = size(S{1},2);
        B = [A{2:sk:end}];
        C = [A{3:sk:end}];
        C2=[];
        for i=1:20
            tmp= Tools.confidence_intervals( [A{sk*(i-1)+1}].' , 95);
            C2(:,i)= tmp(:,2);
        end

        mu_ = mean(B,2).';
        sigm = mean(C2,2).';
        f_ = S{1}{4};
        f_(1)=[];mu_(1)=[];sigm(1)=[];

        fi=(f_>=SYS.analysis_info.f_low & f_<=SYS.analysis_info.f_high);
        mu_oct = Tools.octaveBandMean(mu_(fi),f_(fi),1/6,1000);
        mu_sigm = Tools.octaveBandMean(sigm(fi),f_(fi),1/6,1000);
        SUPP(N_,td_) = mean(mu_oct);
        ConfInt_Low(N_,td_) = mean(mu_sigm);

    end
end

BlockLen = N_vec/1e-3/SYS.signal_info.Fs;
ConfInt_Up = ConfInt_Low;
