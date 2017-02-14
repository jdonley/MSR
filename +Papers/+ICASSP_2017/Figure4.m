%%

SYS = Current_Systems.ICASSP2017_System_A;

SYS.analysis_info.f_high = SYS.signal_info.f_high;

%%
plResult=[];
for a = 1:2
    SYS.signal_info.Nfft = 12 * 1e-3 * SYS.signal_info.Fs;
    if a==1
        SYS.signal_info.time_delay = []; %If empty the time delay will based on the frame length
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
%     Falias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c;
    Falias = 2000;
    
    %%
    if a==1
        col = [1 0 0];
    elseif a==2
        col = [0 0 1];
    end
    yLims = [-20 5];
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
    disp([num2str(mean([D{1:3:end}])), ...
        'dB +-',num2str(mean([D{2:3:end}])), 'dB']);
    mu_ = mean(B,2).';
    sigm = mean(C2,2).';
    f_ = S{1}{4};
    f_(1)=[];mu_(1)=[];sigm(1)=[];
    fi=(f_>=SYS.analysis_info.f_low & f_<=SYS.analysis_info.f_high);
    fsub=f_(fi);
    figure(100);
    ha=area([fsub; fsub].'/1e3, [mu_(fi)+sigm(fi); -2*sigm(fi)].','LineStyle','none','showbaseline','off'); hold on;
    ha(1).FaceColor=col;drawnow;
    ha(2).FaceColor=col;drawnow;
    ha(1).FaceAlpha=0.0;
    ha(2).FaceAlpha=0.3;
    if a==1
        li = '.-';
    elseif a==2
        li = '.--';
    end
    plResult(a) = plot(fsub/1e3,mu_(fi),li,'color',col); hold on;
    set(gca,'XScale','log');
    grid on;grid minor;
    xlim([fix(fsub(1)) ceil(fsub(end))]/1e3);
    ylim(yLims);
    ax=gca;
    ax.XTick=[fix(fsub(1)), 1000 ceil(fsub(end))]/1e3;
    xlabel('Frequency (kHz)');
    ylabel('Suppression (dB)');
    if SYS.signal_info.time_delay == 0
        pred = 'Known Signal';
    elseif ~isempty(SYS.signal_info.time_delay)
        pred = [num2str(SYS.signal_info.time_delay*1e3) 'ms AR'];
    else
        pred = [num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms AR'];
    end
    % title(['Block Size: ' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms   Time Delay: ' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms   Prediction: ' pred])
    title(['Suppression with 12ms Block'])
    set(gca,'FontName','times','FontSize',9);
    ax.YAxis.MinorTick = 'on';
    ax.YAxis.MinorTickValues = yLims(1):yLims(end);
    ax.YAxis.MinorTickValues(1:5:end)=[];
    grid on;
    ax.YMinorGrid = 'on';
    plot([1 1]*Falias/1e3,yLims,'--k');
    hold on;
end

hold on;
ha=area(repmat([Falias,ceil(fsub(end))],2,1).'/1e3, repmat([yLims(1) diff(yLims)],2,1),'LineStyle','none','showbaseline','off'); hold on;
ha(1).FaceColor=[1 1 1];drawnow;
ha(2).FaceColor=[0 0 0];drawnow;
ha(1).FaceAlpha=0.0;
ha(2).FaceAlpha=0.2;
hold off;

legend(plResult,{'Known Signal','Predicted'},'Location','southeast','FontName','times','FontSize',9)

% % % print([SYS.publication_info.DocumentPath '\BS_' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms TD_' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms P_' pred],'-dpng')