%%
for i=flip(1:2)
    res_path={};
    if i==1
        SYS = Current_Systems.IEEETransactions_System_C;
        SYS.signal_info.methods_list((1:end)~=2) = [];
        SYS.signal_info.L_noise_mask(2:end) = [];
        res_path{1} = Hardware_Control.getRealRecordingsPath(SYS);
    elseif i==2
        SYS = Current_Systems.IEEETransactions_System_C;
        SYS.signal_info.methods_list((1:end)~=1) = [];
        SYS.signal_info.L_noise_mask = -inf;
        res_path{1} = Results.getRecordingsPath(SYS);
        SYS = Current_Systems.IEEETransactions_System_C;
        SYS.Main_Setup = SYS.Masker_Setup(1);
        SYS.signal_info.methods_list((1:end)~=2) = [];
        SYS.signal_info.L_noise_mask(2:end) = [];
        res_path{2} = Results.getRecordingsPath(SYS);
    end
    
    B_=[];
    Q_=[];
    for s = 1:numel(res_path)
        files = Tools.getAllFiles( res_path{s} );
        filesB = files(find(~cellfun(@isempty,strfind( files, 'Bright'))));
        filesQ = files(find(~cellfun(@isempty,strfind( files, 'Quiet'))));
        
        B = [];
        Q = [];
        for f = 1:numel(filesB)
            SB = load(filesB{f});
            SQ = load(filesQ{f});
            if i==1
                SB.Rec_Sigs_B = SB.Rec_Sigs_B.';
                SQ.Rec_Sigs_Q = SQ.Rec_Sigs_Q.';
            end
            B = [B; SB.Rec_Sigs_B.'];
            Q = [Q; SQ.Rec_Sigs_Q.'];
        end
        if size(B_,1)~=0
            B(size(B_,1)+1:end,:)=[];
            Q(size(Q_,1)+1:end,:)=[];
        end
        B_(:,:,s) = B;
        Q_(:,:,s) = Q;
    end
    B = sum(B_,3);
    Q = sum(Q_,3);
    %     B = B_(:,:,1);
    %     Q = Q_(:,:,1);
    
    F = 3;
    if i==1;
        b=[];q=[];
        for r = 1:size(B,2)
            b(:,r) = decimate(double(B(:,r)),F);
            q(:,r) = decimate(double(Q(:,r)),F);
        end
    elseif i==2
        b = B;
        q = Q;
    end
    fs_ = 16000;
    Nfft = 2048;
    [Sb,f] = pwelch(b,hamming(Nfft),Nfft/2,Nfft,fs_,'power');
    [Sq,f] = pwelch(q,hamming(Nfft),Nfft/2,Nfft,fs_,'power');
    
    if i==1
        c='r';lw=2;
    elseif i==2
        c='b';lw=2;
        %         Sq = Sq(:,18);
        %         Sb = Sb(:,11);
    end
    plot(f,mean(pow2db(Sb),2) - ...
          mean(pow2db(Sq),2) ...
        ,c,'linewidth',lw);hold on;
    set(gca,'xscale','log');grid off;grid on;grid minor;
    xlim([80 8000]);
    ylim([-10 50])
end
plot([1 1]*250,[-100 100],'k')
plot([1 1]*150,[-100 100],'k')
plot([1 1]*1500,[-100 100],'k')
hold off;

%%
% L=24;
% calibRecPath = 'Z:\+Calibration_Data\+Recordings\Recording_2016-11-23_22.28.wav';
% y = [];
% Rec_Sigs_B = [];
% [ytmp,fs] = audioread(calibRecPath);
% y = [y; ytmp];
% D = reshape(D,[numel(D)/L L]);
% F = 3; d=[];
% for l = 1:L
%     d(:,l) = decimate(D(:,l),F);
% end
% fs_ = fs/F;
% [b,a] = cheby1(6,1,[100 7000]/(fs_/2));
% d__ = filter(b,a,d);
% Nfft = 1024;Sd=[];
% for l = 1:L
% [Sd(:,:,l),frqs,t] = spectrogram(d__(:,l),hamming(Nfft),Nfft/2,Nfft,16000,'power');
% end
% imagesc(frqs,t,((pow2db(mean(abs(Sd(:,:,:)),3)))).');caxis([-5 25])
%
%
% %%
% close all
%  sortMatrix = {'Female_SX306_', ...
% 'Female_SA2_' ...
% 'Male_SX229_' ...
% 'Female1_SX374_' ...
% 'Male_SI1669_' ...
% 'Female1_SX104_' ...
% 'Female1_SA2_' ...
% 'Male_SA2_' ...
% 'Male_SI2299_' ...
% 'Female1_SI1544_' ...
% 'Female_SX216_' ...
% 'Male1_SX5_' ...
% 'Male_SX49_' ...
% 'Male1_SI545_' ...
% 'Female1_SX14_' ...
% 'Male1_SX185_' ...
% 'Female1_SX194_' ...
% 'Female_SX126_' ...
% 'Male_SX319_' ...
% 'Male_SX139_'};
% [~,I] = sort(sortMatrix);
% Iorder(I)=1:numel(sortMatrix);
%
%     SYS = Current_Systems.loadCurrentSRsystem;
%     SYS.signal_info.methods_list(1:end-1) = [];
%     SYS.signal_info.L_noise_mask(2:end) = [];
%     res_path = Hardware_Control.getRealRecordingsPath(SYS);
%
%     calibRecPath = 'Z:\+Calibration_Data\+Recordings\Recording_2016-11-23_22.28.wav';
%
%     files = Tools.getAllFiles( res_path );
%     filesOrig = files(find(~cellfun(@isempty,strfind( files, 'Original'))));
%     filesRec = files(find(~cellfun(@isempty,strfind( files, 'Bright'))));
%
%     filesOrig = sort(filesOrig);filesOrig = filesOrig(Iorder);
%     filesRec= sort(filesRec);filesRec = filesRec(Iorder);
%
%     figure(1);E=[];P=[];
%     for f = 1:numel(filesOrig)
%         y = [];
%         Rec_Sigs_B = [];
%         [ytmp,fs] = audioread(filesOrig{f});
%         y = [y; ytmp];
%         S = load(filesRec{f});
%         Rec_Sigs_B = [Rec_Sigs_B; S.Rec_Sigs_B];
%
%         D = y/rms(y);
%         R = Rec_Sigs_B/rms(Rec_Sigs_B);
%         F = 3;
%         d = decimate(D,F);
%         r = decimate(double(R),F);
%         fs_ = fs/F;
%         [d_,r_]=alignsignals(d,r);
%         [b,a] = cheby1(6,1,[100 7000]/(fs_/2));
%         d__ = filter(b,a,d_);
%         r__ = filter(b,a,r_);
%         d__(numel(r__)+1:end)=[];
%         Nfft = 256*3;
%         [Sd,frqs,t] = spectrogram(d__,hamming(Nfft),Nfft/2,Nfft,fs_,'power');
%         [Sr,frqs,t] = spectrogram(r__,hamming(Nfft),Nfft/2,Nfft,fs_,'power');
%         Id=reshape(1:numel(filesOrig),5,4).';subplot(4,5,Id(f));
%         P(f) = Tools.pesq_mex_vec(d__,r__,fs_);
%         disp(['PESQ: ' num2str(P(f))])
%         Z=(db2mag(-mag2db(abs(Sd.')) + mag2db(abs(Sr.'))));
%         imagesc(frqs,t,imgaussfilt(Z,0.8));
%         title(['PESQ: ' num2str(P(f))])
%         view(2);colormap gray;caxis([2 10])
%         ylim([0.26 3.3])
%         xlim([100 7000])
%         E = [E;Z(numel(t(t<0.26)):end-numel(t(t<0.26)),frqs>150 & frqs<7000)];
%         % figure;pwelch(d__,hamming(1024),512,1024,16000);set(gca,'XScale','log');grid on;grid minor;
%         %     fprintf('\n')
%
%     end
%
%     %%
%     figure(100);
%     imagesc(frqs,t,mag2db(abs(Sd.')));
%     colormap gray;caxis([-10 50]); xlabel('Frequency (Hz)');ylabel('Time (s)');
%     titleTXT = {'Original';['PESQ: ' num2str(P(end)) '   Nfft: ' num2str(Nfft/fs_*1e3) 'ms']};
%     title(titleTXT);
%     print([strrep(strrep([titleTXT{:}],':',''),' ','_') '.png'],'-dpng','-r900');
%
%     fh1=figure(101);
%     imagesc(frqs,t,mag2db(abs(Sr.')));
%     colormap gray;caxis([-10 50]); xlabel('Frequency (Hz)');ylabel('Time (s)');
%     titleTXT = {'Degraded';['PESQ: ' num2str(P(end)) '   Nfft: ' num2str(Nfft/fs_*1e3) 'ms']};
%     title(titleTXT);
%     print([strrep(strrep([titleTXT{:}],':',''),' ','_') '.png'],'-dpng','-r900');
%
% %%
% % SYS = Current_Systems.loadCurrentSRsystem;
% % res_path = Results.getRecordingsPath(SYS);
% %
% % files = Tools.getAllFiles( res_path );
% % files = files(find(~cellfun(@isempty,strfind( files, 'Quiet'))));
% %
% % speech = [];
% % for f = 1:numel(files)
% %     S = load(files{f});
% %     speech = [speech; S.Rec_Sigs_Q.'];
% % end
% %
% % SPECT_LTASS=[];
% % for r = 1:size(speech,2)
% %    [SPECT_LTASS(:,r),frqs] = Tools.LTASS(speech(:,r),SYS.analysis_info.Nfft,SYS.signal_info.Fs);
% % end
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % % Shape the noise to the speech being used
% % [~, SS, f_SS] ...
% %     = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilt_SS( ...
% %     [], SYS.signal_info );
% %
% % % Shape noise spectrum to match quiet zone leakage spectrum
% % [ ~, ~, QZS, QZS_low, f_QZS ] ...
% %     = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_QZS( ...
% %     [], SYS.Masker_Setup(1), SYS.Main_Setup, SYS.system_info, SYS.signal_info );
% %
% % % Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
% % cheby_order = 9; Rp = 1;
% % f_cutoff = Broadband_Tools.getAliasingFrequency(SYS.Masker_Setup(1)) * SYS.signal_info.c / (2*pi);
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % if any(frqs - f_SS), error('Frequency vectors don''t match.');end
% %
% % f_l = SYS.analysis_info.f_low;
% % f_c = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c;
% % f_h = SYS.analysis_info.f_high;
% %
% % f_band = (f_l<=frqs & frqs<=f_c);
% % f_band2 = (f_l<=f_QZS & f_QZS<=f_c);
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % P=[];
% % for r = 1:size(speech,2)
% %     [P(:,r),f_]  = Tools.octaveBandMean(SPECT_LTASS(f_band,r),frqs(f_band),1/6);
% % end
% % octSpace = 1/6;
% % [sp,f1] = Tools.octaveBandMean(SS(f_band),frqs(f_band),octSpace);
% % [q1,f2] = Tools.octaveBandMean(QZS(f_band2),f_QZS(f_band2),octSpace);
% %  q2     = Tools.octaveBandMean(QZS_low(f_band2),f_QZS(f_band2),octSpace);
% %  lp     = 1 ./ sqrt( 1 + (10^(Rp/10)-1)*chebyshevT(cheby_order,(f1/f_cutoff)).^2);
% % if any(f1 - f2), error('Frequency vectors don''t match.');end
% % f1([1 end])=[];
% % sp([1 end])=[];
% % q1([1 end])=[];
% % q2([1 end])=[];
% % lp([1 end])=[];
% % P([1 end],:)=[];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % wh = ones(size(f1)) * mean(P(:));
% % p = ones(size(f1)) * mean(P(:)) ./ (f1/mean(f1));
% % spq2lp = sp.*q2.*lp;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % H = [sp;q2;lp;wh;p;spq2lp]; E=[]; Eps=[];
% % Eps_min = ones(1,size(H,1))*1e4;
% % for h = 1:size(H,1)
% %     H(h,:) = Tools.Correlated_Normalisation(mean(P,2),H(h,:));
% %
% %     for g = db2mag(-30:.1:30)
% %         [E(:,h),Eps(h)]=Tools.COSHdist(H(h,:)*g,P);
% %         Eps_min(h) = min([Eps(h),Eps_min(h)]);
% %     end
% % end
% %
% % Eps_min
% % %%
% % % f=400;
% % % f2 = 480;
% % % T = 1/f;
% % % fs=16000;
% % % t = 0:1/fs:0.032;t(end)=[];
% % % x = cos(2*pi*f*t + pi/4) + 0.5*cos(2*pi*f2*t + pi/3);% + rand(1,length(t))*0.2;
% % % w = hanning(numel(x)).';
% % % x_ = x;%.*w;
% % % f_ = linspace(0,fs/2,numel(x_)/2+1);f_(end)=[];
% % %
% % % td = 0.015;
% % %
% % % X = fft(x_);
% % % XX = X(1:end/2) .* exp(1i*2*pi*f_*td);
% % % y = ifft( [XX, 0, conj(XX(end:-1:2))] );
% % %
% % % x__=x_;
% % % for ii=1:numel(x__)/2
% % %     a = lpc(x__);
% % %     xn = -sum( a(2:end) .* x__(end:-1:end-numel(a)+2) );
% % %     x__ = [x__, xn];
% % % end
% % %
% % % f__ = linspace(0,fs/2,numel(x__)/2+1);f__(end)=[];
% % % X = fft(x__);
% % % XX = X(1:end/2) .* exp(1i*2*pi*f__*td);
% % % y2 = ifft( [XX, 0, conj(XX(end:-1:2))] );
% % %
% % % plot(x__);hold on;
% % % plot(x_);hold on;
% % % % plot(y,'k');
% % % % plot(y2,'r');
% % % hold off;
% %
% % %%
% % % load('Z:\_RECOR~1\_LINEA~1.565\_1GENE~1\_0BX_0~1\_VSRC_~1.581\_DATAB~2\_10X10~1\32REC_~1\_150HZ~1\ESS_100Hz-7kHz_10sec_1secTail_1secHead_Bright.mat')
% % % T = Rec_Sigs_B;
% % % load('Z:\_RECOR~1\_2LINE~1\_36GEN~1.074\_0BX_0~1\_VSRC_~1.581\_DATAB~2\_10X10~1\32REC_~1\_150HZ~1\ESS_100Hz-7kHz_10sec_1secTail_1secHead_Bright.mat')
% % % M = Rec_Sigs_B;
% % % load('M:\MSR\+Miscellaneous\+TestAudio_Files_InvFilts\inverseESS.mat');
% % %
% % %
% % % figure(1);
% % % for x = 1:8
% % %     for y = 1:4
% % %         p=(y-1)*8+x;
% % %         IRt(p,:) = Tools.extractIR( T(p,:), invY);
% % %         IRm(p,:) = Tools.extractIR( -M(p,:), invY);
% % %
% % %         subplot(4,8,p);
% % %         plot(IRt(p,:)); hold on;
% % %         plot(IRm(p,:)); hold off;
% % %         [~,I]=max([IRt(p,:)]);
% % %         xlim([-1 1]*10 + I);
% % %         ylim([-1 1]*max([IRt(p,:),IRm(p,:)]));
% % %     end
% % % end
% %
% % %%
% % SYS = Current_Systems.loadCurrentSRsystem;
% % N_vec = (4:4:32)*1e-3*SYS.signal_info.Fs;
% % td_vec = {0,[]};
% % mu_octM = zeros(numel(N_vec),numel(td_vec));
% % mu_sigmM = mu_octM;
% %
% % for td_ = 1:numel(td_vec)
% %     for N_ = 1:numel(N_vec)
% %         SYS.signal_info.Nfft = N_vec(N_);
% %         SYS.signal_info.time_delay = td_vec{td_};
% %         SYS.signal_info.zeropadtime = SYS.signal_info.Nfft / SYS.signal_info.Fs;
% %         SYS.system_info.LUT_frequencies = (SYS.signal_info.Nfft + SYS.signal_info.zeropadtime * SYS.signal_info.Fs)/2;
% %         SYS.system_info.LUT_resolution = [num2str(SYS.system_info.LUT_frequencies) 'f' ...
% %                               SYS.system_info.sc ...
% %                               num2str(SYS.system_info.LUT_weights) 'w'];
% %         SYS.Main_Setup = SYS.Main_Setup(1);
% %         SYS.signal_info.method = SYS.signal_info.methods_list{end};
% %         load([Results.getResultsPath(SYS) 'Suppression_Results.mat']);
% %         A = [S{:}]; sk = size(S{1},2);
% %         B = [A{2:sk:end}];
% %         C = [A{3:sk:end}];
% %         C2=[];
% %         for i=1:20
% %             tmp= Tools.confidence_intervals( [A{sk*(i-1)+1}].' , 95);
% %             C2(:,i)= tmp(:,2);
% %         end
% %
% %         mu_ = mean(B,2).';
% %         sigm = mean(C2,2).';
% %         f_ = S{1}{4};
% %         f_(1)=[];mu_(1)=[];sigm(1)=[];
% %
% %         fi=(f_>=SYS.analysis_info.f_low & f_<=SYS.analysis_info.f_high);
% %         mu_oct = Tools.octaveBandMean(mu_(fi),f_(fi),1/6,1000);
% %         mu_sigm = Tools.octaveBandMean(sigm(fi),f_(fi),1/6,1000);
% %         mu_octM(N_,td_) = mean(mu_oct);
% %         mu_sigmM(N_,td_) = mean(mu_sigm);
% %
% %     end
% % end
% %
% % figure(101);
% % M = repmat(N_vec/1e-3/SYS.signal_info.Fs,2,1).';
% % errorbar(M,mu_octM,mu_sigmM,mu_sigmM,'-^','linewidth',1)
% % ax = gca;
% % ax.Children(1).Color = [0.9 0 0];
% % ax.Children(1).Marker = 'x';
% % ax.Children(2).Color = [0.1 0.3 0.8];
% % ax.XTickMode = 'manual';
% % ax.YTickMode = 'manual';
% % ax.ZTickMode = 'manual';
% % ax.XMinorTickMode = 'manual';
% % ax.YMinorTickMode = 'manual';
% % ax.ZMinorTickMode = 'manual';
% % ax.XLimMode = 'manual';
% % ax.YLimMode = 'manual';
% % ax.ZLimMode = 'manual';
% %      ax.XAxis.MinorTick = 'on';
% %      ax.XAxis.MinorTickValues = 0:M(end,1)+10;
% % %      ax.XAxis.MinorTickValues(1:4:33)=[];
% %      ax.YAxis.MinorTick = 'on';
% %      ax.YAxis.MinorTickValues = -30:0;
% % %      ax.YAxis.MinorTickValues(1:5:31)=[];
% %      grid on;
% %      ax.XMinorGrid = 'on';
% %      ax.YMinorGrid = 'on';
% % xlim(M([1 end])+[-2 2])
% % ylim([-20 0]+[-2 0])
% % set(gca,'xtick',M(:,1))
% % set(gca,'ytick',-20:5:0)
% % xlabel('Block Length (ms)')
% % ylabel('Suppression (dB)')
% % title('Synthesis and Prediction Trade-off')
% % legend({'Known Signal','Predicted'},'Location','southwest','FontName','times','FontSize',9)
% % set(gca,'FontName','times','FontSize',9);
% % set(gcf, 'PaperPositionMode', 'auto');
% %
% %
% %
% % %%
% % SYS = Current_Systems.loadCurrentSRsystem;
% % SYS.Main_Setup = SYS.Main_Setup(1);
% % SYS.signal_info.method = SYS.signal_info.methods_list{end};
% % load([Results.getResultsPath(SYS) 'Suppression_Results.mat'])
% % A = [S{:}]; sk = size(S{1},2);
% % B = [A{2:sk:end}];
% % C = [A{3:sk:end}];
% % C2=[];
% % for i=1:20
% %     tmp= Tools.confidence_intervals( [A{sk*(i-1)+1}].' , 95);
% %     C2(:,i)= tmp(:,2);
% % end
% % D = [A{sk:sk:end}];
% % disp([num2str(mean([D{1:3:end}])), ...
% %     'dB +-',num2str(mean([D{2:3:end}])), 'dB']);
% % mu_ = mean(B,2).';
% % sigm = mean(C2,2).';
% % f_ = S{1}{4};
% % f_(1)=[];mu_(1)=[];sigm(1)=[];
% % fi=(f_>=SYS.analysis_info.f_low & f_<=SYS.analysis_info.f_high);
% % fsub=f_(fi);
% % figure(100);
% % ha=area([fsub; fsub].'/1e3, [mu_(fi)+sigm(fi); -2*sigm(fi)].','LineStyle','none','showbaseline','off'); hold on;
% % ha(1).FaceColor=[1 1 1]*1;drawnow;
% % ha(2).FaceColor=[1 1 1]*0;drawnow;
% % ha(1).FaceAlpha=0.0;
% % ha(2).FaceAlpha=0.3;
% % plot(fsub/1e3,mu_(fi),'.-','color',[1 1 1]*0); hold on;
% % set(gca,'XScale','log');
% % grid on;grid minor;
% % xlim([fix(fsub(1)) 2000]/1e3);
% %  ylim([-15 0]);
% % ax=gca;
% % ax.XTick=[fix(fsub(1)), 1000 2000]/1e3;
% % xlabel('Frequency (kHz)');
% % ylabel('Suppression (dB)');
% % if SYS.signal_info.time_delay == 0
% %     pred = 'Known Signal';
% % elseif ~isempty(SYS.signal_info.time_delay)
% %     pred = [num2str(SYS.signal_info.time_delay*1e3) 'ms AR'];
% % else
% %     pred = [num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms AR'];
% % end
% % % title(['Block Size: ' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms   Time Delay: ' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms   Prediction: ' pred])
% % title(['Suppression with Optimal Predicted Block'])
% % set(gca,'FontName','times','FontSize',9);
% %      ax.YAxis.MinorTick = 'on';
% %      ax.YAxis.MinorTickValues = -15:0;
% %      ax.YAxis.MinorTickValues(1:5:end)=[];
% %      grid on;
% %      ax.YMinorGrid = 'on';
% %
% % hold off;
% % % print([SYS.publication_info.DocumentPath '\BS_' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms TD_' num2str(SYS.signal_info.Nfft/SYS.signal_info.Fs*1e3) 'ms P_' pred],'-dpng')
% % %%
% % % %
% % % f=1500;
% % % c=343;
% % % fmid=1000;
% % % d=c/fmid/12;
% % %
% % % nT=5;
% % % x = linspace(-c/f*nT,c/f*nT,500);
% % % x2 = x+d;
% % % [xx,yy]=meshgrid(x);
% % % xx2=meshgrid(x2,x);
% % % [tt,rr]=cart2pol(xx,yy);
% % % [tt2,rr2]=cart2pol(xx2,yy);
% % %
% % % % r = abs(x);
% % % % r2 = abs(x2);
% % %
% % % k = 2*pi*f/c;
% % % a1 = besselh(0,k*rr);
% % %
% % % % compshift = exp(1i*pi*(1/2-2*d*f/c));
% % % % b1 = besselh(0,k*rr) * compshift;
% % % % c1 = besselh(0,k*rr2)*exp(1i*(1+2*d*f/c)*pi) * compshift;
% % %
% % % b1 = besselh(0,k*rr) * exp( -1i*pi/2) * exp( 1i*pi*2*d*f/c);
% % % c1 = besselh(0,k*rr2) * exp( -1i*pi/2) * exp( 1i*pi);
% % %
% % % figure(11)
% % % surf( real(a1) ,'linestyle','none');view(2)
% % % figure(1)
% % % surf( real(b1+c1) ,'linestyle','none');view(2)
% % %
% % % figure(2)
% % % plot(x,real(a1(end/2,:)),'k'); hold on;
% % % plot(x,real(-a1(end/2,:)),':k'); hold on;
% % % plot(x,real(b1(end/2,:)),'--r');
% % % plot(x,real(c1(end/2,:)),'--m');
% % % plot(x,real(b1(end/2,:)+c1(end/2,:)),'b');
% % % grid on; grid minor;
% % % hold off
% % % axis([-c/f*nT c/f*nT -2 2])
% % % %%
% % % % clc;
% % % % SpkrSigs=audioread('Z:\+Speaker_Signals\+CIRCLEarray_1.3mPerpDist_180degCentre\+24Genelec8010ASpkrs_4.2616mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_2w\+150Hz-8000Hz_-InfdB_100weight__method_NoMask\ESS_100Hz-7kHz_10sec_1secTail_1secHead_24Ch.WAV');
% % % % Orig=audioread('Z:\+Speaker_Signals\+CIRCLEarray_1.3mPerpDist_180degCentre\+24Genelec8010ASpkrs_4.2616mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_2w\+150Hz-8000Hz_-InfdB_100weight__method_NoMask\ESS_100Hz-7kHz_10sec_1secTail_1secHead_Original.WAV');
% % % % load('M:\MSR\+Miscellaneous\+TestAudio_Files_InvFilts\inverseESS.mat');
% % % %
% % % % load('Z:\+Recordings\+CIRCLEarray_1.3mPerpDist_180degCentre\+24Genelec8010ASpkrs_4.2616mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_2w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_100weight__method_NoMask\ESS_100Hz-7kHz_10sec_1secTail_1secHead_Bright.mat')
% % % %
% % % % IR = Tools.extractIR( SpkrSigs(:,10), invY);
% % % % IRo = Tools.extractIR( Orig, invY);
% % % %
% % % % IRr=[];
% % % % for i = 1:32
% % % % irtmp = Tools.extractIR( Rec_Sigs_B(i,:), invY);
% % % % [~,IRr(i,:)] = alignsignals(IRo,irtmp,8000,'truncate');
% % % % end
% % % % IRr = IRr.';
% % % %
% % % % [v,I]=findpeaks(IR/max(IR),'NPeaks',2,'MinPeakHeight',0.01,'MinPeakDistance',500)
% % % % figure(1);
% % % % subplot(1,3,1);plot(IR/max(IR));
% % % % subplot(1,3,2);plot(IRr/max(IRr(:)));
% % % % subplot(1,3,3);plot(IRo/max(IRo));
% % % % diff(I)
% % % %
% % % %
% % % % %%
% % % % % width=length(obj.Soundfield_d);
% % % % % x = ((1:width) - width/2) / obj.res;
% % % % % y = x;
% % % % % [xx,yy] = meshgrid(x,y);
% % % % % X = complex(xx,yy);
% % % % %
% % % % % obj.Soundfield_d = weight .* besselh(0, k * abs(X .* exp(1i.*(pi-Phi_src)) + R_src));
% % % % %
% % % % %
% % % % % %%
% % % % % A = fft(Tools.fconv(data,flip(data)));
% % % % % for r = 1:32
% % % % % B(:,r) = fft(Tools.fconv(data,flip(Rec_Sigs_B(r,:).')));
% % % % % end
% % % % %
% % % % % figure(1)
% % % % % F = linspace(0,8000,numel(A)/2);
% % % % % plot(F,unwrap(angle(A(1:end/2)))); hold on;
% % % % % plot(F,mean(unwrap(angle(B(1:end/2,:))),2)); hold on;
% % % % % set(gca,'XScale','log'); grid on; grid minor;xlim([0.1 8]*1e3);
% % % % % hold off
% % % % %
% % % % % figure(4)
% % % % % plot(F,mag2db(abs(A(1:end/2)))); hold on;
% % % % % plot(F,mean(mag2db(abs(B(1:end/2,:))),2)+40); hold on;
% % % % % set(gca,'XScale','log'); grid on; grid minor;xlim([0.1 8]*1e3);ylim([90 110]);
% % % % % hold off
% % % % % %
% % % % % %
% % % % % % figure(2)
% % % % % % plot(linspace(0,8000,numel(B)/2),unwrap(angle(A(1:end/2)))-unwrap(angle(B(1:end/2)))); hold on;
% % % % % % hold off;
% % % % %
% % % % % figure(3)
% % % % % [d,R]=alignsignals(data,Rec_Sigs_B(1,:));
% % % % % plot(d./rms(d))
% % % % % hold on
% % % % % plot(R./rms(R))
% % % % % hold off
% % % % %
% % % % %
% % % % % %%
% % % % %
% % % % % f = 100:10:8000;
% % % % % K = f*2*pi/343;
% % % % %
% % % % %
% % % % % a=0;b=0;
% % % % % MM = nan(numel(K),2*ceil(K(end)*1.0)+1);
% % % % % for k_ = 1:numel(K)
% % % % %     k=K(k_);
% % % % %     M = ceil(k*1.0);
% % % % %     Mv = -M:M;
% % % % %     for m_ = 1:numel(Mv)
% % % % %         m = Mv(m_);
% % % % %         a(k_,m_) =  1./besselh(m, k)  ;
% % % % %
% % % % %         b(k_,m_) =  1./besselh(m, 1000*2*pi/343 )  ;
% % % % %
% % % % %         MM(k_,m_) = m;
% % % % %     end
% % % % %     disp(k_);
% % % % % end
% % % % % a = sum(a,2);
% % % % % b = sum(b,2);
% % % % %
% % % % % b_ = unwrap(angle( b ));
% % % % % b_norm = (b_ - b_(f==1000));
% % % % % plot(f, b_norm ); hold on
% % % % % plot(f, unwrap(angle( a )) ); hold on
% % % % %
% % % % % KK = repmat(K.',1,2*ceil(K(end)*1.0)+1);
% % % % % T = 1./(besselh(MM, KK ));
% % % % % T(isnan(T))=0; T=1./sum(T,2);
% % % % % plot(f, unwrap(angle( a.* T )));
% % % % % set(gca,'XScale','log'); grid on; grid minor;
% % % % % hold off
% % % % %
% % % % % %%
% % % % % % e=[];
% % % % % % for f = 100:8000
% % % % % % k=f*2*pi/343;
% % % % % % M = ceil(k);
% % % % % % e(f-99)=0;
% % % % % % for m = -M:M
% % % % % % e(f-99) = e(f-99) + abs( sqrt(abs(besselh(0, k))) );
% % % % % % end
% % % % % % end
% % % % %
% % % % % figure(1)
% % % % % [p,f]=pwelch(Rec_Sigs_B.',hamming(1024),512,1024,16000,'power');
% % % % % [po,f_]=pwelch(data,hamming(1024),512,1024,16000,'power');
% % % % % p2 = abs(fft(Rec_Sigs_B.'));
% % % % % po2 = abs(fft(data));
% % % % % %plot(f_,mag2db( mean(db2mag(pow2db(p)),2) )+39.5/2,'k');hold on;
% % % % % plot(linspace(0,8000,size(p2,1)/2),mean(mag2db( p2(1:end/2,:).' ),1)+40,'r');hold on;
% % % % % plot(linspace(0,8000,size(po2,1)/2),mag2db( po2(1:end/2) ),'m');hold on;
% % % % % %plot(f_,mag2db( mean(db2mag(pow2db(po)),2) ),'b');hold on;
% % % % %
% % % % % % plot(100:8000,mag2db(e));hold on;
% % % % % % plot(100:8000,mag2db( abs(besselh(0, (100:8000)*2*pi/343))  ));hold on;
% % % % %
% % % % %
% % % % % % plot(f_,mag2db( mean(db2mag(pow2db(p)),2) .* abs(besselh(0, (f_)*2*pi/343)) ),'k');hold on;
% % % % %
% % % % % hold off;
% % % % %
% % % % % set(gca,'XScale','log');
% % % % % xlim([0.1 8]*1e3);
% % % % % grid on; grid minor;
% % % % %
% % % % %
% % % % %
% % % % % figure(2)
% % % % % hold off
% % % % % for i = 1:32
% % % % % [d,R]=alignsignals(data,Rec_Sigs_B(i,:));
% % % % % a2 = angle(fft(R));
% % % % % ao2 = angle(fft(d));
% % % % % % plot(unwrap(ao2)); hold on;
% % % % % % plot(unwrap(a2(:,1)));
% % % % % hold on;
% % % % % plot(linspace(0,8000,floor(numel(ao2)/2)),unwrap(ao2(1:end/2)) );
% % % % % plot(linspace(0,8000,floor(numel(ao2)/2)),unwrap(a2(1:floor(numel(ao2)/2))) );
% % % % % end
% % % % % xlim([0.1 8]*1e3);
% % % % % grid on; grid minor;
% % % % % hold off
% % % % %
% % % % % figure(3);
% % % % % [d,R]=alignsignals(data,Rec_Sigs_B(1,:));
% % % % % plot(d./rms(d)); hold on
% % % % % plot(R./rms(R)); hold on
% % % % %
% % % % % %%
% % % % % % close all;
% % % % % % f1 = 20; %Hz
% % % % % % f2 = 20000; %Hz
% % % % % % fs = 48000; %Hz
% % % % % % c = 343; %metres/sec
% % % % % % len = 30; %sec
% % % % % % tail = 5; %sec
% % % % % %
% % % % % % y = Tools.synthSweep(len+tail,fs,f1,f2,tail*fs);
% % % % % % invFFT = Tools.invSweepFFT(y,f1,f2,fs);
% % % % % % y = [zeros(tail*fs,1); y(:)];
% % % % % % %audiowrite('sinesweep_20Hz-20kHz_30sec_5secTail.WAV',y,fs);
% % % % % %
% % % % % % %[y2, fs2] = audioread('result.WAV');
% % % % % % hold on;
% % % % % % delay_avg = []; eq=[];
% % % % % %
% % % % % % for dist = 1:2
% % % % % %     for s = 1:2
% % % % % %         for e_ = 1
% % % % % % %  for spkr = 2
% % % % % % %  for rec = 1
% % % % % %
% % % % % % %     data = audioread(['+Results\+Test_Recordings\+23_11_2015\sinesweep_20Hz-20kHz_30sec_5secTail_spkr' num2str(spkr) '_dist' num2str(dist) '_rec' num2str(rec) '.wav']);
% % % % % %     switch s
% % % % % %         case 1
% % % % % %             Shift = 'no';
% % % % % %         case 2
% % % % % %             Shift = 'with';
% % % % % %     end
% % % % % %     EQ = 'no';
% % % % % %     data = audioread(['+Results\+Test_Recordings\+24_11_2015\distance' num2str(dist) '_' EQ 'EQ_' Shift 'Shift.wav']);
% % % % % %
% % % % % %     y2 = data(1:(length(y)-(tail-1)*fs));
% % % % % %
% % % % % % %     f1_=f1;
% % % % % % %     X = fft(y2);
% % % % % % %     % Frequencies
% % % % % % % M = length(X);
% % % % % % % NumPts = M/2 + 1;
% % % % % % % freqs = linspace(0, fs/2, NumPts);
% % % % % % %
% % % % % % %     W = FrequencyResponse((freqs>=f1_) & (freqs<=f2));
% % % % % % %      I=find(freqs>1000,1,'first');
% % % % % % %      W =  (W ./ W(I));
% % % % % % %      %W = smooth(W,1000);
% % % % % % %
% % % % % % %      W = -mag2db(W);
% % % % % % %
% % % % % % % %     plot(freqs((freqs>=f1_) & (freqs<=f2)),W);
% % % % % % % %     grid on;
% % % % % % % %     set(gca,'YScale','lin');
% % % % % % % %     set(gca,'XScale','log');
% % % % % % %
% % % % % % % % Weight Levels
% % % % % % % W = [linspace(0, W(1), length(freqs(freqs<f1_))), ...
% % % % % % %      W(:)', ...
% % % % % % %      linspace(W(end), 0, length(freqs(freqs>f2)))];
% % % % % % %
% % % % % % %
% % % % % % %
% % % % % % %
% % % % % % % % Apply magnitude weighting
% % % % % % % X(1:NumPts) = X(1:NumPts) .* db2mag(W(:));
% % % % % % %
% % % % % % % % Apply conjugation for negative frequency side of spectrum
% % % % % % % X(NumPts+1:M) = conj(X(M/2:-1:2));
% % % % % % %
% % % % % % % % Time Domain Transform
% % % % % % % data_ = ifft(X); % Inverse Fast Fourier Transform
% % % % % % %
% % % % % % % % prepare output vector y
% % % % % % % data_ = real(data_(1:M));
% % % % % % %
% % % % % % % % remove DC
% % % % % % % data_ = data_(:) - mean(data_);
% % % % % % %     y2=data_;
% % % % % %
% % % % % %     %y2 = y2 ./ max(abs(y2(:)));
% % % % % %
% % % % % %     [ir,irinv] = Tools.extractIR(y2,invFFT);
% % % % % %
% % % % % %     ir_ = ir((tail*fs):end);
% % % % % %
% % % % % %     irfull = [irinv; ir];
% % % % % %     Zr = fft(irfull);
% % % % % %     fr = abs(Zr( 1:((end/2)-1) ));
% % % % % %
% % % % % %     %eq((spkr-1)*5+rec,:) = fr;
% % % % % %
% % % % % %     f = linspace(0,fs/2,length(Zr)/2);f(end)=[];
% % % % % %
% % % % % %     figure(1); hold on;
% % % % % %     plot(f,(mag2db(fr)));
% % % % % %     set(gca,'XScale','log');
% % % % % %     xlim([20,20000]);
% % % % % %     ylim([-80,0]);
% % % % % %     xlabel('Frequency (Hz)');
% % % % % %     ylabel('Magnitude (dB)');
% % % % % %     grid on;
% % % % % %
% % % % % %     figure(2); hold on;
% % % % % %     plot((1:length(ir_))/fs*c, ir_);
% % % % % %     [pkv,pk]=max(ir_(:));
% % % % % %     plot(pk/fs*c,pkv,'or');
% % % % % %     xlim([0, 3.0]);
% % % % % %     xlabel('Distance (metres)');
% % % % % %     grid on; grid minor;
% % % % % %     delay_avg(dist,s) = pk;
% % % % % % %  end
% % % % % % %  end
% % % % % % end
% % % % % % end
% % % % % % end
% % % % % % hold off;
% % % % % % delay_ = diff(mean(delay_avg,2));
% % % % % % distance_ = delay_/fs * c; % metres
% % % % % % figure(2);
% % % % % % title(['Distance moved: ' num2str(round(distance_*100,2)) ' centimetres']);
% % % % % %
% % % % % % %%
% % % % % % % L = 24;
% % % % % % % b = setup.Loudspeaker_Locations(L,2);
% % % % % % % c = 0.6;
% % % % % % % A = 270-setup.Loudspeaker_Locations(L,1)/pi*180;
% % % % % % % a = sqrt(b^2 + c^2 - 2*b*c*cosd(A));
% % % % % % %
% % % % % % % B = asind( sind(A) * b / a );
% % % % % % %
% % % % % % % maxAngle = 90-B
% % % % % %
% % % % % %
% % % % % % % c = 343;
% % % % % % % phiL = 90/180*pi;
% % % % % % % R = 0.9;
% % % % % % %
% % % % % % % L = 1:24;
% % % % % % %
% % % % % % % f = (L-1)*c / (2*R*phiL);
% % % % % % %
% % % % % % % plot(f,L)
% % % % % % %
% % % % % % % f = f(1):f(end);
% % % % % % %
% % % % % % % L = ceil( ((phiL/pi*180) / 360) * 2 * ( ceil( (f/c *2*pi) * R )) + 1 );
% % % % % % %
% % % % % % % hold on;
% % % % % % %
% % % % % % % plot(f,L)
% % % % % % %
% % % % % % % hold off
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % % %%
% % % % % % % load('LUT_Weight_vs_Frequency_15deg_512f_32w.mat')
% % % % % % % load('LUT_Weight_vs_Frequency_0deg_512f_32w_zones_swapped.mat')
% % % % % % % load('LUT_Weight_vs_Frequency_15deg_512f_32w_compact.mat')
% % % % % % % load('LUT_Weight_vs_Frequency_0deg_512f_32w_zones_swapped_compact.mat')
% % % % % % % load('Database_90degB_-90degQ_0deg_512f_32w.mat')
% % % % % % % load('Database_-90degB_90degQ_-60deg_512f_32w.mat')
% % % % % %
% % % % % % %LUT_MagDiff = (Bright_Sample__Weight_Vs_Frequency ./ Quiet_Sample__Weight_Vs_Frequency);
% % % % % % % LUT_MagDiff = mag2db(Acoustic_Contrast__Weight_Vs_Frequency);
% % % % % % %
% % % % % % % Number_of_Weights = 256;
% % % % % % % len = 4096;
% % % % % % % Fs=16000;
% % % % % % % noise_freqs = linspace(0, Fs/2, len/2 + 1);
% % % % % % % noise_freqs = noise_freqs(noise_freqs>=min(Frequencies) & noise_freqs<=max(Frequencies));
% % % % % % % weights = [0 logspace(log10( 1e-2), log10(  1e4), Number_of_Weights - 1) ];
% % % % % % %
% % % % % % % LUT_MagDiff_interp = interp2(Frequencies,Weights,LUT_MagDiff,noise_freqs',weights,'spline');
% % % % % % %
% % % % % % % [val,I]=max(LUT_MagDiff_interp);
% % % % % % % noise_weights = weights(I);
% % % % % % %
% % % % % % %
% % % % % % % figure(1);
% % % % % % % surf(noise_freqs,weights, LUT_MagDiff_interp,'linestyle','none');
% % % % % % % view(2);
% % % % % % % set(gca,'XScale','log');
% % % % % % % set(gca,'YScale','log');
% % % % % % % set(gca,'TickDir','out');
% % % % % % % ylabel('Quiet Zone Weight');
% % % % % % % xlabel('Frequency (Hz)');
% % % % % % % c=colorbar;
% % % % % % % ylabel(c,'Zone Contrast (dB)')
% % % % % % % %caxis([0 1]);
% % % % % % % hold on;
% % % % % % %
% % % % % % %  plot3(noise_freqs, noise_weights, ones(1,length(noise_freqs)).*max(abs(val(:)))*1.1 ,'k');
% % % % % % %  hold off;
% % % % % % %
% % % % % % %
% % % % % % % W = -mag2db(Quiet_Sample__Weight_Vs_Frequency);
% % % % % % % W_ = permute( Tools.interpVal_2D(W, Frequencies, Weights, noise_freqs, noise_weights, 'spline'), [2 1]);
% % % % % % % W_=W_(:);
% % % % % % %
% % % % % % % %hold on
% % % % % % % figure(2);
% % % % % % % plot(noise_freqs, val/2);
% % % % % % % set(gca,'XScale','log');
% % % % % % % set(gca,'YScale','lin');
% % % % % % % set(gca,'TickDir','out');
% % % % % % % ylabel('Mag Difference');
% % % % % % % xlabel('Frequency (Hz)');
% % % % % % % %ylim([0 1]);
% % % % % % % grid on;
% % % % % % %
% % % % % % % hold on
% % % % % % % %figure(3);
% % % % % % % plot(noise_freqs, W_);
% % % % % % % set(gca,'XScale','log');
% % % % % % % set(gca,'YScale','lin');
% % % % % % % set(gca,'TickDir','out');
% % % % % % % ylabel('Level (dB)');
% % % % % % % xlabel('Frequency (Hz)');
% % % % % % % grid on;