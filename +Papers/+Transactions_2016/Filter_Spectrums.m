clc;clear;close all;

%%
SYS = Current_Systems.IEEETransactions_System_F;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Signal_Length = 180 * SYS.signal_info.Fs;

WGN = v_addnoise( zeros(Signal_Length,1), SYS.signal_info.Fs, -Inf); % White noise

% Shape the noise to the speech being used
[SSN, SS, f_SS] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilt_SS( ...
    WGN, SYS.signal_info );

% Quiet zone leakage spectrum
[ ~, ~, QZS, QZS_low, f_QZS ] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_ZoneSpect( 'Quiet', ...
    [], SYS.Masker_Setup(1), SYS.Main_Setup, SYS.system_info, SYS.signal_info );

% Bright zone second order leakage spectrum
[ ~, ~, BZS, BZS_low, f_BZS ] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_ZoneSpect( 'Quiet', ...
    [], SYS.Masker_Setup(1), SYS.Masker_Setup(1), SYS.system_info, SYS.signal_info );

if f_QZS ~= f_BZS, error('Frequency vectors don''t match'); end
ZS_ratio = sqrt( QZS ./ BZS );
ZS_ratio_low = sqrt( QZS_low ./ BZS_low );
SSN_ZS = Tools.ArbitraryOctaveFilt(SSN, ZS_ratio, f_QZS, SYS.signal_info.Nfft, SYS.signal_info.Fs, SYS.signal_info.OctaveBandSpace);
SSN_ZS_low = Tools.ArbitraryOctaveFilt(SSN, ZS_ratio_low, f_QZS, SYS.signal_info.Nfft, SYS.signal_info.Fs, SYS.signal_info.OctaveBandSpace);

% Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
cheby_order = 9; Rp = 1;
SSN_ZS_LP = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_AliasCtrl( ...
    SSN_ZS_low, SYS.Masker_Setup(1), SYS.signal_info, cheby_order, Rp );

Noise_Signal1 = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_PassbandNorm( ...
    SSN_ZS, SSN_ZS_LP, SYS.Masker_Setup(1), SYS.signal_info ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WGN = v_addnoise( zeros(Signal_Length,1), SYS.signal_info.Fs, -Inf); % White noise

% Shape the noise to the speech being used
[SSN, SS, f_SS] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilt_SS( ...
    WGN, SYS.signal_info );


% Shape noise spectrum to match quiet zone leakage spectrum
[ SSN_QZS, SSN_QZS_low, QZS, QZS_low, f_QZS ] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_ZoneSpect( 'Quiet', ...
    SSN, SYS.Masker_Setup(1), SYS.Main_Setup, SYS.system_info, SYS.signal_info );

% Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
cheby_order = 9; Rp = 1;
SSN_QZS_LP = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_AliasCtrl( ...
    SSN_QZS_low, SYS.Masker_Setup(1), SYS.signal_info, cheby_order, Rp );

% Adjust the power of the signal in the passband to match the non-filtered version (power normalisation (equalisation))
Noise_Signal2 = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_PassbandNorm( ...
    SSN_QZS, SSN_QZS_LP, SYS.Masker_Setup(1), SYS.signal_info ); 

%
% [A,fA] = pwelch(SSN_ZS,rectwin(SYS.analysis_info.Nfft),0,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');A=sqrt(A);
% [B,fB] = pwelch(SSN_ZS_LP,rectwin(SYS.analysis_info.Nfft),0,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');B=sqrt(B);
% plot([fA,fB],mag2db([A,B]));set(gca,'XScale','log'); grid on; grid minor;xlim([0.10 8]*1e3);ylim([-90 -40]-0);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_cutoff = Broadband_Tools.getAliasingFrequency(SYS.Masker_Setup(1))/2/pi*SYS.signal_info.c;

spect_sp = SS;
frqs_sp = f_SS;
spect = QZS;
spect_lowpass = QZS_low;
frqs = f_QZS;
Input_Signal_filt_I = SSN_QZS_LP;
Input_Signal_filt_IB = SSN_ZS_LP;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
win_=rectwin(SYS.analysis_info.Nfft);
ovlap = 0;

res_path = Results.getRecordingsPath(SYS);
files = Tools.getAllFiles( res_path );
files = Tools.keepFilesFromFolder(files,SYS.signal_info.speech_filepath);
filesQ = files(find(~cellfun(@isempty,strfind( files, 'Quiet'))));
filesB = files(find(~cellfun(@isempty,strfind( files, 'Bright'))));
speech = [];speechB = [];
for f = 1:numel(filesQ)
    S = load(filesQ{f});
    speech = [speech; S.Rec_Sigs_Q.'];
end
for f = 1:numel(filesB)
    S = load(filesB{f});
    speechB = [speechB; S.Rec_Sigs_B.'];
end
SPECT_LTASS=[]; SPECT_LTASS_B=[];
for r = 1:size(speech,2)
   [SPECT_LTASS(:,r),freqs] = pwelch(speech(:,r),win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');SPECT_LTASS(:,r)=sqrt(SPECT_LTASS(:,r));
   [SPECT_LTASS_B(:,r),freqsB] = pwelch(speechB(:,r),win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');SPECT_LTASS_B(:,r)=sqrt(SPECT_LTASS_B(:,r));
end

%%
SYS.signal_info.L_noise_mask = 0;
SYS.signal_info.method = SYS.signal_info.methods_list{SYS.signal_info.methods_list_masker(end)};
mask_path = Results.getRecordingsPath(SYS.Masker_Setup(end),SYS.system_info.LUT_resolution, SYS.Room_Setup, SYS.signal_info, SYS.system_info.Drive);
files = Tools.getAllFiles( mask_path );
files = files(find(~cellfun(@isempty,strfind( files, 'Quiet'))));
masker = [];
for f = 1:numel(files)
    S = load(files{f});
    masker = [masker; S.Rec_Sigs_Q.'];
end
MASKER_LTASS=[];
for r = 1:size(masker,2)
   [MASKER_LTASS(:,r),freqs] = pwelch(masker(:,r),win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');MASKER_LTASS(:,r)=sqrt(MASKER_LTASS(:,r));
end

%%
octSpace = SYS.signal_info.OctaveBandSpace;
P=[];P_B=[];M=[];
for r = 1:size(speech,2)
    [P(:,r),f_]   = Tools.octaveBandMean( SPECT_LTASS(:,r),freqs,octSpace);
    [P_B(:,r),f_] = Tools.octaveBandMean( SPECT_LTASS_B(:,r),freqsB,octSpace);
    [M(:,r),f_m]  = Tools.octaveBandMean(MASKER_LTASS(:,r),freqs,octSpace);
end
%
Srec = mean(P,2).*db2mag(0);
SrecB = mean(P_B,2).*db2mag(0);
Srecmask = mean(M,2).*db2mag(-15);

[Ufinal_I,fp] = pwelch(Input_Signal_filt_I,win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');Ufinal_I=sqrt(Ufinal_I);
[Ufinal_IB,fp] = pwelch(Input_Signal_filt_IB,win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');Ufinal_IB=sqrt(Ufinal_IB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Figure Output Settings
DocumentPath = SYS.publication_info.DocumentPath;
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = 88.9/10;% + 6.35/10 + 88.9/10; %IEEE full text width
aspect_ratio = 4/5;
FontSize = 8;
Font = 'Times';
lineWid = 0.5;
LegendHeightScaleFactor = 1.1;

% Latex Fonts
FF = SYS.publication_info.LaTeX_FontFamily; % Set fonts for latex interpreter
FFnums = SYS.publication_info.LaTeX_NumbersFontFamily; % Set number fonts for latex interpreter
FS = SYS.publication_info.FontSize;
latexFontSettings = ['{\fontfamily{' FF '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];
latexNumFontSettings = ['{\fontfamily{' FFnums '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];

% Start Plotting
ha = tightPlots(3,1,plot_width,[2.0 0.9], [0.1 0.1], 1, 1,'centimeters');

axLblXpos = 150;
xlims = [0.1 8];
ylims = [-50 0];
numXTicks = 3;
xticks = [xlims(1) 1.0 xlims(end)];
XTlbls = mat2cell(num2str(xticks'),ones(size(num2str(xticks'),1),1),size(num2str(xticks'),2));
XTlbls = cellfun(@(x) strrep(x,' ',''), XTlbls, 'un', 0);
XTlbls = cellfun(@(x) [latexNumFontSettings x '}'],XTlbls,'un',0);
figParams = {'XScale','log', ...
    'TickLabelInterpreter', 'latex', ...
    'XLim', xlims, ...
    'XTick', xticks, ...
    'XTickLabelMode', 'manual', ...
    'XTickLabel', XTlbls, ...
    'YTickMode', 'manual', ...
    'TickDir', 'both', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'XMinorGrid','on', ...
    'YMinorGrid','off', ...
    'FontName',Font, ...
    'FontSize',FontSize};
axes(ha(1));

[S1,f1]=Tools.octaveBandMean(spect_sp,frqs_sp,octSpace);
[S2,f2]=Tools.octaveBandMean(spect,frqs,octSpace);
S3=Tools.octaveBandMean(spect_lowpass,frqs,octSpace);
S4=Tools.octaveBandMean(BZS,frqs,octSpace);
S5=Tools.octaveBandMean(BZS_low,frqs,octSpace);
S6 = 1 ./ sqrt( 1 + (10^(Rp/10)-1)*chebyshevT(cheby_order,(f1/f_cutoff)).^2);

ms = 4;
f_c = f_cutoff*[1 1]/1e3;
pl_Ku1= plot(f_c,ylims,'--k', 'LineWidth',lineWid,'markersize',ms); hold on;
pl_S1 = plot(f1/1e3,mag2db(S1),'-',  'Color',[0 0 0],       'LineWidth',lineWid,'markersize',ms); hold on;
pl_S2 = plot(f2/1e3,mag2db(S2),'--', 'Color',[0.8 0.5 0.5], 'LineWidth',lineWid,'markersize',ms); hold on;
pl_S3 = plot(f2/1e3,mag2db(S3),'--', 'Color',[1 0 0],       'LineWidth',lineWid,'markersize',ms); hold on;
pl_S4 = plot(f2/1e3,mag2db(S4),'-.', 'Color',[0.5 0.5 0.8],'LineWidth',lineWid,'markersize',ms); hold on;
pl_S5 = plot(f2/1e3,mag2db(S5),'-.', 'Color',[0 0 1],      'LineWidth',lineWid,'markersize',ms); hold on;
pl_S6 = plot(f1/1e3,mag2db(S6),':',  'Color',[0 0 0],        'LineWidth',lineWid,'markersize',ms); hold off;

text(axLblXpos/1e3,ylims(end)-7,'(A)','interpreter','latex','horizontalalignment','center');
set(gca,figParams{:},'YLim',ylims);% axis square
set(gca,'XTickLabel',{});
yticks = [ylims(1):10:ylims(end)].';
Tlbls= [repmat(latexNumFontSettings,size(yticks,1),1) ...
    num2str(yticks) ...
    repmat('}',size(yticks,1),1)];
[m,n]=size(Tlbls);
YTickLbls=mat2cell(Tlbls,ones(m,1),n);
set(gca,'YTick',yticks);
set(gca,'YTickLabel',YTickLbls);
set(gca,'fontname','fixedwidth');
title('Filter Spectra','fontweight','normal','fontname','times', 'FontSize',FontSize);
% xlabel('Frequency (kHz)','fontname','times', 'FontSize',FontSize);
% ylabel('Magnitude (dB)','fontname','times', 'FontSize',FontSize);
[leg1, legIcons] = legend([pl_S1;pl_S2;pl_S3;pl_S4;pl_S5;pl_S6],...
    {'$H^{(\mathrm{sp})}$';...
    '$H^{(\mathrm{q})}$';...
    '$H^{(\mathrm{q''})}$';...
    '$H^{(\mathrm{b})}$';...
    '$H^{(\mathrm{b''})}$';...
    '$H^{(\mathrm{lp})}$';},...
    'box','off','color','none', ...
    'interpreter','latex','location','eastoutside','fontsize',FontSize);
leglines = findobj(legIcons,'Type','line');
legShift = -min(diff(reshape([leglines(1:2:length(leglines)).XData],2,[])).*2/5);
for ll = 1:2:length(leglines)
    leglines(ll).XData(1) = leglines(ll).XData(1) - legShift;
end
leg1.Units = 'centimeters';
NewHeight = leg1.Position(4) * LegendHeightScaleFactor;
leg1.Position(2) = leg1.Position(2) - (NewHeight - leg1.Position(4))/2;
leg1.Position(4) = NewHeight;




axes(ha(2));

H_sp = S1;
H_q1 = S3([1 end-length(S1)+2:end]);
H_lp = S6;
H_b1 = S5([1 end-length(S1)+2:end]);

f1band = (f1>SYS.signal_info.f_low & f1<f_cutoff);
f_band = (f_>SYS.signal_info.f_low & f_<f_cutoff);
SQuiet = Srec;
SBright = SrecB;

lambda = 0;
SQuietMask = H_sp .* (H_q1.^(1-lambda)) ./ (H_b1.^(lambda)) .* H_lp * db2mag(0);
lambda = 1/2;
SIdealMask = H_sp .* (H_q1.^(1-lambda)) ./ (H_b1.^(lambda)) .* H_lp * db2mag(0);
lambda = 1;
SBrightMask = H_sp .* (H_q1.^(1-lambda)) ./ (H_b1.^(lambda)) .* H_lp * db2mag(0);

SQuietMaskB = SQuietMask .* H_b1;
SIdealMaskB = SIdealMask .* H_b1;
SBrightMaskB = SBrightMask .* H_b1;


pl_0 = plot(f_/1e3, mag2db(SQuiet / mean(SQuiet(f_band)) ) ,'-','color',[0 0 0 1],'LineWidth',lineWid); hold on
pl_2 = plot(f1/1e3, mag2db(SQuietMask / mean(SQuietMask(f1band)) ),'--','color',[1 0 0 1],'LineWidth',lineWid); hold on
pl_3 = plot(f1/1e3, mag2db(SIdealMask / mean(SIdealMask(f1band)) ),'.-','color',[0 0 1 1],'LineWidth',lineWid); hold on
pl_4 = plot(f1/1e3, mag2db(SBrightMask / mean(SBrightMask(f1band)) ),'-.','color',[0 0.5 0 1],'LineWidth',lineWid); hold on

hold off
ylims = ylims+20;
text(axLblXpos/1e3,ylims(end)-7,'(B)','interpreter','latex','horizontalalignment','center');
set(gca,figParams{:},'YLim',ylims);% axis square
set(gca,'XTickLabel',{' '});
yticks = [ylims(1):10:ylims(end)].';
Tlbls= [repmat(latexNumFontSettings,size(yticks,1),1) ...
    num2str(yticks) ...
    repmat('}',size(yticks,1),1)];
[m,n]=size(Tlbls);
YTickLbls=mat2cell(Tlbls,ones(m,1),n);
set(gca,'YTick',yticks);
set(gca,'YTickLabel',YTickLbls);
set(gca,'fontname','fixedwidth');
% title('Filtered Masker','fontweight','normal','fontname','times', 'FontSize',FontSize);
% xlabel('Frequency (kHz)','fontname','times', 'FontSize',FontSize);
ylabel('Magnitude (dB)','fontname','times', 'FontSize',FontSize);
[leg2, legIcons] = legend([pl_0;...
    pl_2;...
    pl_3;...
    pl_4;],...
    {'${\bar{P}}^{(\mathrm{sp,q})}$';...
    ['$H^{(\mathrm{{\mathcal{IB}},lp})},$' 10 '$\,\,\lambda{\grave{}}=0$'];...
    ['$H^{(\mathrm{{\mathcal{IB}},lp})},$' 10 '$\,\,\lambda{\grave{}}=0.5$'];...
    ['$H^{(\mathrm{{\mathcal{IB}},lp})},$' 10 '$\,\,\lambda{\grave{}}=1$'];},...
    'box','off','color','none', ...
    'interpreter','latex','location','eastoutside','fontsize',FontSize);
leglines = findobj(legIcons,'Type','line');
set(leglines,'LineWidth',lineWid); % Force linewidth in legend to all be the same for better display
legShift = -min(diff(reshape([leglines(1:2:length(leglines)).XData],2,[])).*2/5);
for ll = 1:1:length(leglines)
    leglines(ll).XData(1) = leglines(ll).XData(1) - legShift;
    if numel(leglines(ll).XData) == 1
        leglines(ll).XData(1) = leglines(ll).XData(1) + legShift/2;
    end
end
leg2.Units = 'centimeters';
NewHeight = leg2.Position(4) * LegendHeightScaleFactor;
leg2.Position(2) = leg2.Position(2) - (NewHeight - leg2.Position(4))/2;
leg2.Position(4) = NewHeight;






axes(ha(3));

pl_1 = plot(f_/1e3, mag2db(SBright / mean(SBright(f_band)) ) ,'-','color',[0 0 0 1],'LineWidth',lineWid); hold on
pl_5 = plot(f1/1e3, mag2db(SQuietMaskB / mean(SQuietMaskB(f1band)) ),'--','color',[1 0 0 1],'LineWidth',lineWid); hold on
pl_6 = plot(f1/1e3, mag2db(SIdealMaskB / mean(SIdealMaskB(f1band)) ),'.-','color',[0 0 1 1],'LineWidth',lineWid); hold on
pl_7 = plot(f1/1e3, mag2db(SBrightMaskB / mean(SBrightMaskB(f1band)) ),'-.','color',[0 0.5 0 1],'LineWidth',lineWid); hold on


hold off
ylims = ylims;
text(axLblXpos/1e3,ylims(end)-7,'(C)','interpreter','latex','horizontalalignment','center');
set(gca,figParams{:},'YLim',ylims);% axis square
yticks = [ylims(1):10:ylims(end)].';
Tlbls= [repmat(latexNumFontSettings,size(yticks,1),1) ...
    num2str(yticks) ...
    repmat('}',size(yticks,1),1)];
[m,n]=size(Tlbls);
YTickLbls=mat2cell(Tlbls,ones(m,1),n);
set(gca,'YTick',yticks);
set(gca,'YTickLabel',YTickLbls);
set(gca,'fontname','fixedwidth');
% title('Filtered Masker','fontweight','normal','fontname','times', 'FontSize',FontSize);
xlabel('Frequency (kHz)','fontname','times', 'FontSize',FontSize);
% ylabel('Magnitude (dB)','fontname','times', 'FontSize',FontSize);
[leg3, legIcons] = legend([pl_1;...
    pl_5;...
    pl_6;...
    pl_7;],...
    {'${\bar{P}}^{(\mathrm{sp})}$';...
    ['$H^{(\mathrm{{\mathcal{IB}},b,lp})},$' 10 '$\,\,\lambda{\grave{}}=0$'];...
    ['$H^{(\mathrm{{\mathcal{IB}},b,lp})},$' 10 '$\,\,\lambda{\grave{}}=0.5$'];...
    ['$H^{(\mathrm{{\mathcal{IB}},b,lp})},$' 10 '$\,\,\lambda{\grave{}}=1$'];},...
    'box','off','color','none', ...
    'interpreter','latex','location','eastoutside','fontsize',FontSize);
leglines = findobj(legIcons,'Type','line');
set(leglines,'LineWidth',lineWid); % Force linewidth in legend to all be the same for better display
legShift = -min(diff(reshape([leglines(1:2:length(leglines)).XData],2,[])).*2/5);
for ll = 1:1:length(leglines)
    leglines(ll).XData(1) = leglines(ll).XData(1) - legShift;
    if numel(leglines(ll).XData) == 1
        leglines(ll).XData(1) = leglines(ll).XData(1) + legShift/2;
    end
end
leg3.Units = 'centimeters';
legShiftDist = -legShift*leg3.Position(3);

NewHeight = leg3.Position(4) * LegendHeightScaleFactor;
leg3.Position(2) = leg3.Position(2) - (NewHeight - leg3.Position(4))/2;
leg3.Position(4) = NewHeight;








% export figure for publication
h=gcf;
tightfig(h);
h.Units = 'centimeters';
h.Color = 'w';
fig_pos = h.Position;
h.Position = [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio];
h.PaperUnits = 'centimeters';
h.PaperSize = [plot_width plot_width/aspect_ratio];
l1P = leg1.Position;
leg1.Position = l1P + [-0.0 0 0 0];
l2P = leg2.Position;
leg2.Position = [l1P(1) l2P(2) l2P(3) l2P(4)] + [-0.0 0 0 0];
l3P = leg3.Position;
leg3.Position = [l1P(1) l3P(2) l3P(3) l3P(4)] + [-0.0 0 0 0];


set(ha(1),'Units','centimeters');
set(ha(2),'Units','centimeters');
set(ha(3),'Units','centimeters');

ax1 = get(ha(1));ax2 = get(ha(2));ax3 = get(ha(3));
set(ha(1),'OuterPosition',[0,ax1.Position(2),h.Position(3)-leg1.Position(3)+legShiftDist,ax1.Position(4)]);
set(ha(2),'OuterPosition',[0,ax2.Position(2),h.Position(3)-leg1.Position(3)+legShiftDist,ax2.Position(4)]);
set(ha(3),'OuterPosition',[0,ax3.Position(2),h.Position(3)-leg1.Position(3)+legShiftDist,ax3.Position(4)]);

ax1 = get(ha(1));ax2 = get(ha(2));ax3 = get(ha(3));
set(ha(2),'Position',[ax1.Position(1),ax2.Position(2),ax1.Position(3),ax1.Position(4)]);
set(ha(3),'Position',[ax1.Position(1),ax3.Position(2),ax1.Position(3),ax1.Position(4)]);

ax1 = get(ha(1));ax2 = get(ha(2));ax3 = get(ha(3));
leg1.Position(1) = h.Position(3)- leg1.Position(3);
leg2.Position(1) = leg1.Position(1);
leg3.Position(1) = leg1.Position(1);
leg1.Position(2) = ax1.Position(2) + (ax1.Position(4)-leg1.Position(4))/2;
leg2.Position(2) = ax2.Position(2) + (ax2.Position(4)-leg2.Position(4))/2;
leg3.Position(2) = ax3.Position(2) + (ax3.Position(4)-leg3.Position(4))/2;


%
tightfig(h);

% ht=suptitle('Acoustic Contrast for Reproduction Methods');
% ht.FontName = 'times';
% ht.FontSize = 10;
% ht.Position = ht.Position + [0 -0.08 0];

if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
% print(['-d' print_fmt], [DocumentPath '\MSR_PL_HYB_AcousticContrast_matlab.pdf']);
  export_fig([DocumentPath '\Masker_Filters_' num2str(SYS.Main_Setup.Dimensionality) 'D.pdf']);

%
Tools.MiKTeX_FNDB_Refresh;