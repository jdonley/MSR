function [Freqs,SPL,ConfInt_Low,ConfInt_Up] = import_SoundPressureLevels_vs_Freq( SYS )
%import_SoundPressureLevels_vs_Freq Import numeric data from a text file as column vectors.
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up]
%   = import_SoundPressureLevels_vs_Freq(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up]
%   = import_SoundPressureLevels_vs_Freq(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up] = import_Suppression_Reverb('Suppression_Results.mat',1, 180);
%

%%
SPL = [];
ConfInt_Low=[];
ConfInt_Up=[];
%%
rec_types = SYS.signal_info.recording_type;
failCnt = 0;
for rt = 1:numel(SYS.signal_info.recording_type)
    SYS.signal_info.recording_type = rec_types(rt);
        %%
        SYS.signal_info.method = 'Sweep'; % This software framework uses an exponential sine sweep method to accurately compute the SPL per frequency.
        if ~exist([Results.getResultsPath(SYS) 'SPL_Results.mat'],'file')
            failCnt = failCnt+1;
            continue;
        end
        if failCnt == numel(rec_types)
            Freqs = -1;
            SPL = ['''' SYS.signal_info.recording_type{:} ''' results file not found.'];
            return;
        end
        S = load([Results.getResultsPath(SYS) 'SPL_Results.mat']);S=S.S;
        A = [S{:}]; sk = size(S{1},2);
        
        B = [A{2}];
        CB = [A{3}];
        
        Q = [A{5}];
        CQ = [A{6}];
        
        f__ = A{7}; f__(1)=[]; f_(:,rt) = f__;
        fi=(f_(:,rt)>=SYS.analysis_info.f_low & f_(:,rt)<=SYS.signal_info.f_high);
        Falias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c;
        fiMean=(f_(:,rt)>=SYS.signal_info.f_low & f_(:,rt)<=Falias);
        meanB = mean(mean(B(fiMean,:)));
        
        SPLB(:,rt) = mean(B,2).' - meanB;
        SPLQ(:,rt) = mean(Q,2).' - meanB;
        ConfIntB_Low(:,rt) = mean(CB,2).';
        ConfIntQ_Low(:,rt) = mean(CQ,2).';
        
end

SPLB(1,:)=[];SPLQ(1,:)=[];ConfIntB_Low(1,:)=[];ConfIntQ_Low(1,:)=[];
Freqs=f_(fi) / 1e3; %KiloHertz
SPLB(~fi,:)=[];
SPLQ(~fi,:)=[];
ConfIntB_Low(~fi,:)=[];
ConfIntQ_Low(~fi,:)=[];

SPL(:,:,1) = SPLB;
SPL(:,:,2) = SPLQ;
ConfInt_Low(:,:,1) = ConfIntB_Low;
ConfInt_Low(:,:,2) = ConfIntQ_Low;
ConfInt_Up = ConfInt_Low;



