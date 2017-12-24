%%
%  load('Z:\+Soundfield_Database\+3m_SpkrDia\+16Spkrs_180DegArc_original\LUT_Weight_vs_Frequency_15deg_256f_128w.mat')

%%
fh = figure(54321);
close(fh);

plot_width = 88.9/10;% + 6.35/10 + 88.9/10; %IEEE full text width (cm)

fh = figure(54321);
axH = tightPlots(2, 1, ...
    plot_width, ...
    [1 1/3], ...
    [0.5 0.5], ...
    [1 1], ...
    [1 1], ...
    'centimeters');

axes(axH(1)); ax1 = gca; ax1.Units = 'centimeters';
surf(Frequencies/1e3,...
    mag2db(Weights),...
    mag2db(abs(Bright_Sample__Weight_Vs_Frequency)),...
    'lines','no');
view(2);
set(gca,'xscale','log')
colormap(gray);
caxis([-50 -15])
xlim([0.15 8]);
ylim([-40 80]);
ax1.TickDir = 'both';
ax1.XTick = [0.15 1 8];
ax1.XTickLabel = {};
cbH1 = colorbar;
cbH1.Units = ax1.Units;
cbH1.Location = 'manual';
cbH1.Position = [0.25 + sum(ax1.Position([1 3])), ...
                ax1.Position(2), ...
                0.5, ...
                ax1.Position(4)];

axes(axH(2)); ax2 = gca; ax2.Units = 'centimeters';
surf(Frequencies/1e3,...
    mag2db(Weights),...
    mag2db(abs(Quiet_Sample__Weight_Vs_Frequency)),...
    'lines','no');
view(2);
set(gca,'xscale','log')

colormap(gray);
caxis([-50 -15]);
xlim([0.15 8]);
ylim([-40 80]);
ax2.TickDir = 'both';
ax2.XTick = [0.15 1 8];
cbH2 = colorbar;
cbH2.Units = ax2.Units;
cbH2.Location = 'manual';
cbH2.Position = [0.25 + sum(ax2.Position([1 3])), ...
                ax2.Position(2), ...
                0.5, ...
                ax2.Position(4)];
cbH2.Label.String = 'SPL (dB)';
cbH2.Label.Units = ax2.Units;

xlabel('Frequency (kHz)');
ylabel('Quiet Zone Weight (dB)');
ax2.YLabel.Units = ax2.Units;


cbH2.Label.Position(2) = [sum(cbH1.Position([2 4])) - cbH2.Position(2)]/2
ax2.YLabel.Position(2) = [sum(ax1.Position([2 4])) - ax2.Position(2)]/2


tightfigadv;

