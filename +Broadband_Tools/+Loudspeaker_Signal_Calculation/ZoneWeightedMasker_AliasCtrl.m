function ZoneWeightedMasker_AliasCtrl( Signal_Length, SYS_or_LUT_resolution, Noise_Mask_dB, weight, setup, multizone_setup )
%ZONEWEIGHTEDMASKER_ALIASCTRL Generates spatially weighted noise masker
%signals with a pre-processed noise signal.

SYS_type = 'Current_Systems.SR_System';

%%
plotPublicationFig = false; % true to plot publication figure
figN = 1111; % figure number
if ishandle(figN)
    figure(figN+1)
    return;
end

%% Setup Variables
if nargin == 2
    if isa(SYS_or_LUT_resolution,SYS_type)
        SYS = SYS_or_LUT_resolution;        
        signal_info = SYS.signal_info;
        system_info = SYS.system_info;
        setup = SYS.Masker_Setup;
        multizone_setup = SYS.Main_Setup;
        if any(size(signal_info.L_noise_mask)~=1)
           warning('Different number of noise levels may not be supported.'); 
        end
    else
        error(['Second input argument must be of type: ' SYS_type]);
    end
elseif nargin == 6
    system_info.Drive = 'Z:\'; % Database drive (storage drive)
    system_info.LUT_resolution = SYS_or_LUT_resolution;
    
    signal_info.c = 343; % Speed of sound in metres/sec
    signal_info.Fs = 16000; % Sampling frequency
    signal_info.Nfft = 1024;% Number of fft components
    signal_info.overlap = 0.5;
    signal_info.f_low  = 150;  % Hz
    signal_info.f_high = 8000; % Hz
    signal_info.L_noise_mask = Noise_Mask_dB; % dB
    signal_info.weight = weight;
    signal_info.method = 'ZoneWeightMaskerAliasCtrl';
end

signal_info.input_filename = 'Masker';

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, system_info.LUT_resolution, system_info.Drive, 'new');


%% Generate Input Signal
WGN = v_addnoise( zeros(Signal_Length,1), signal_info.Fs, -Inf); % White noise

%% Shape the noise to the speech being used
[SSN, SS, f_SS] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilt_SS( ...
    WGN, signal_info );

%% Shape noise spectrum to match quiet zone leakage spectrum
[ SSN_QZS, SSN_QZS_low, QZS, QZS_low, f_QZS ] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_QZS( ...
    SSN, setup, multizone_setup, system_info, signal_info );

%% Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
cheby_order = 9; Rp = 1;
SSN_QZS_LP = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_AliasCtrl( ...
    SSN_QZS_low, setup, signal_info, cheby_order, Rp );

%% Adjust the power of the signal in the passband to match the non-filtered version (power normalisation (equalisation))
Noise_Signal = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_PassbandNorm( ...
    SSN_QZS, SSN_QZS_LP, setup, signal_info, -35 ); % Hope that -35dB RMS level is low enough to avoid clipping upon saving

%% Plot Filter Responses for Publication
if plotPublicationFig
    figure(figN);
    f_cutoff = Broadband_Tools.getAliasingFrequency(setup) * signal_info.c / (2*pi);
    plotPublicationFigure(SS,f_SS,QZS,QZS_low,f_QZS,Rp,cheby_order,f_cutoff,SSN_QZS_LP);
    figure(figN+1)
    return
end

%% Compute the loudspeaker signals for the additive zone weighted noise
Loudspeaker_Signals = Broadband_Tools.Loudspeaker_Signal_Calculation.getMSRLoudspeakerSignals( ...
    Noise_Signal, setup, signal_info, system_info )...
    *  db2mag(signal_info.L_noise_mask);


%% Once we have the speaker signals we should save them for later use as .wav files
Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    [], [], [], [], ...
    signal_info.Fs );


end



function plotPublicationFigure(spect_sp,frqs_sp,spect,spect_lowpass,frqs,Rp,cheby_order,f_cutoff,Input_Signal_filt)
% Figure Output Settings
DocumentPath = 'tex\latex\TRANS2016';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = 88.9/10;% + 6.35/10 + 88.9/10; %IEEE full text width
aspect_ratio = 6/3;
FontSize = 9;
Font = 'Times';
lineWid = 1.0;

% Start Plotting
ha = tightPlots(1,2,plot_width,[1 1.5], [0 1.5], 1, 1,'centimeters');
xlims = [0.1 10];
numXTicks = 3;
xticks = logspace(log10(xlims(1)),log10(xlims(end)),numXTicks);
XTlbls = mat2cell(num2str(xticks'),ones(size(num2str(xticks'),1),1),size(num2str(xticks'),2));
XTlbls = cellfun(@(y) strrep(y,' ',''), XTlbls, 'UniformOutput', false);
figParams = {'XScale','log', ...
    'XLim', xlims, ...
    'XTick', xticks, ...
    'XTickLabelMode', 'manual', ...
    'XTickLabel', XTlbls, ...
    'XGrid','on', ...
    'YGrid','on', ...
    'fontname',Font, ...
    'FontSize',FontSize};
axes(ha(1));
[S1,f1]=Tools.octaveBandMean(spect_sp,frqs_sp,1/6);
pl_S1 = plot(f1/1e3,mag2db(S1),'k','LineWidth',lineWid);
set(gca,figParams{:},'YLim',[-60 0]-30);axis square
hold on;
[S2,f2]=Tools.octaveBandMean(spect,frqs,1/6);
S3=Tools.octaveBandMean(spect_lowpass,frqs,1/6);
pl_S2 = plot(f2/1e3,mag2db(S2),'-','Color',[0.5 0.5 0.5],'LineWidth',lineWid); hold on;
pl_S3 = plot(f2/1e3,mag2db(S3),'--k','LineWidth',lineWid); hold on;
hold off
set(gca,figParams{:},'YLim',[-60 0]); axis square
hold on;
S4 = 1 ./ sqrt( 1 + (10^(Rp/10)-1)*chebyshevT(cheby_order,(f1/f_cutoff)).^2);
pl_S4 = plot(f1/1e3,mag2db(S4),'-.k','LineWidth',lineWid);
set(gca,figParams{:},'YLim',[-90 0]); axis square
set(gca,'fontname','fixedwidth');
title('Masker Filter Spectrums','fontweight','normal','fontname','times', 'FontSize',FontSize);
xlabel('Frequency (kHz)','fontname','times', 'FontSize',FontSize);ylabel('Magnitude (dB)','fontname','times', 'FontSize',FontSize);
[leg1, legIcons] = legend([pl_S1;pl_S2;pl_S3;pl_S4],...
    {'$H^{(sp)}$';'$H^{(q)}$';'$H^{\prime(q)}$';'$H^{(lp)}$';},...
    'box','on','color','none', ...
    'interpreter','latex','location','southwest','fontsize',8);
leglines = findobj(legIcons,'Type','line');
legShift = -min(diff(reshape([leglines(1:2:length(leglines)).XData],2,[])).*1/5);
for ll = 1:2:length(leglines)
    leglines(ll).XData(2) = leglines(ll).XData(2) + legShift;
end
leg1.Units = 'points';
[leg1,legIcons] = gridLegend([pl_S1;pl_S2;pl_S3;pl_S4],2,...
    'existinglegend',leg1,'existinglegendicons',legIcons,...
    'legendwidth',leg1.Position(3).*(1+legShift).*0.9,...
    'location','southoutside'); leg1.Units = 'centimeters';
hold off;
axes(ha(2));
Sfinal = S1 .* S3([1 end-length(S1)+2:end]) .* S4;
[Ufinal,fp] = pwelch(Input_Signal_filt);
pl_Ufinal = plot(fp/pi*8, pow2db(Ufinal),'color',[0.5 0.5 1],'LineWidth',lineWid); hold on
pl_Sfinal = plot(f1/1e3, mag2db(Sfinal),'k','LineWidth',lineWid); hold off
set(gca,figParams{:},'YLim',[-90 0]-50); axis square
set(gca,'fontname','fixedwidth');
title('Filtered Masker','fontweight','normal','fontname','times', 'FontSize',FontSize);
xlabel('Frequency (kHz)','fontname','times', 'FontSize',FontSize);ylabel('Magnitude (dB)','fontname','times', 'FontSize',FontSize);
[leg2, legIcons] = legend([pl_Ufinal;pl_Sfinal;],...
    {'$\tilde{U}_{\mathrm{a}}^{(sp,q,lp)}$';'$H^{(sp,q,lp)}$';},...
    'box','on','color','none', ...
    'interpreter','latex','location','southwest','fontsize',8);
leglines = findobj(legIcons,'Type','line');
legShift = -min(diff(reshape([leglines(1:2:length(leglines)).XData],2,[])).*1/5);
for ll = 1:2:length(leglines)
    leglines(ll).XData(2) = leglines(ll).XData(2) + legShift;
end
leg2.Units = 'points';
[leg2,legIcons] = gridLegend([pl_Ufinal;pl_Sfinal;],1,...
    'existinglegend',leg2,'existinglegendicons',legIcons,...
    'legendwidth',leg2.Position(3).*(1+legShift).*1.0,...
    'location','southoutside'); leg2.Units = 'centimeters';
%     figure(2);
%     hold on
%     pwelParams = {[], ...
%         0, ...
%         [], ...
%         signal_info.Fs, ...
%         'power'};
%     pwelch(Input_Signal_SSN); hold on
%     pwelch(Input_Signal_SSN_QZshaped); hold on
%     pwelch(Input_Signal_SSN_QZshaped_lowpass); hold on
%     pwelch(Input_Signal); hold on
%     hold off
%     ax = gca;
%     pwel = get(ax,'Children');
%     arrayfun(@(i) set(pwel(i),'Color',ax.ColorOrder(mod(i-1,7)+1,:)),1:length(pwel))
%     set(gca,figParams{:},'XLim',[0.1 10]/10,'YLim',[-90 0]-60);

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
leg1.Position = l1P + [0 0 0 0];
l2P = leg2.Position;
leg2.Position = [l2P(1) l1P(2) l2P(3) l2P(4)] + [0 -0.05 0 0];
tightfig(h);

% ht=suptitle('Acoustic Contrast for Reproduction Methods');
% ht.FontName = 'times';
% ht.FontSize = 10;
% ht.Position = ht.Position + [0 -0.08 0];

if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
% print(['-d' print_fmt], [DocumentPath '\MSR_PL_HYB_AcousticContrast_matlab.pdf']);
export_fig([DocumentPath '\Masker_Filters.pdf']);

end