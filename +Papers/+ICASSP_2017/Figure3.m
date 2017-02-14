%%

SYS = Current_Systems.ICASSP2017_System_A;

% c = SYS.signal_info.c;
% Falias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup(1))/2/pi*343;
% SYS.analysis_info.f_high = Falias;

%%
N_vec = (4:4:32)*1e-3*SYS.signal_info.Fs;
td_vec = {0,[]};
mu_octM = zeros(numel(N_vec),numel(td_vec));
mu_sigmM = mu_octM;

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
        mu_octM(N_,td_) = mean(mu_oct);
        mu_sigmM(N_,td_) = mean(mu_sigm);

    end
end

figure(101);
M = repmat(N_vec/1e-3/SYS.signal_info.Fs,2,1).';
errorbar(M,mu_octM,mu_sigmM,mu_sigmM,'-^','linewidth',1)
ax = gca;
ax.Children(1).Color = [0.9 0 0];
ax.Children(1).Marker = 'x';
ax.Children(2).Color = [0.1 0.3 0.8];
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XMinorTickMode = 'manual';
ax.YMinorTickMode = 'manual';
ax.ZMinorTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues = 0:M(end,1)+10;
%      ax.XAxis.MinorTickValues(1:4:33)=[];
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = -30:0;
%      ax.YAxis.MinorTickValues(1:5:31)=[];
     grid on;
     ax.XMinorGrid = 'on';
     ax.YMinorGrid = 'on';
xlim(M([1 end])+[-2 2])
ylim([-20 0]+[-4 0])
set(gca,'xtick',M(:,1))
set(gca,'ytick',-20:5:0)
xlabel('Block Length (ms)')
ylabel('Suppression (dB)')
title('Synthesis and Prediction Trade-off')
legend({'Known Signal','Predicted'},'Location','southwest','FontName','times','FontSize',9)
set(gca,'FontName','times','FontSize',9);
set(gcf, 'PaperPositionMode', 'auto');