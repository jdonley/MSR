%%
% load('Z:\+Soundfield_Database\+3m_SpkrDia\+16Spkrs_180DegArc_original\LUT_Weight_vs_Frequency_15deg_256f_128w.mat')

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

axes(axH(1)); ax = gca; ax.Units = 'centimeters';
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
ax.TickDir = 'both';
ax.XTick = [0.15 1 8];
ax.XTickLabel = {};
cbH = colorbar;
cbH.Units = ax.Units;
cbH.Location = 'manual';
cbH.Position = [0.25 + sum(ax.Position([1 3])), ...
                ax.Position(2), ...
                0.5, ...
                ax.Position(4)];


axes(axH(2)); ax = gca;
surf(Frequencies/1e3,...
    mag2db(Weights),...
    mag2db(abs(Quiet_Sample__Weight_Vs_Frequency)),...
    'lines','no');
view(2);
set(gca,'xscale','log')

colormap(gray);
colorbar('location','eastoutside');
caxis([-50 -15])
xlim([0.15 8]);
ylim([-40 80]);
ax.TickDir = 'both';
ax.XTick = [0.15 1 8];

xlabel('Frequency (kHz)');
ylabel('Quiet Zone Weight (dB)');


tightfigadv;

