clc;
clear;
close all;
% N Planewave effect on LUT Comparison

LUT_resolution = '512f_256w';
angle_pw = 15;
loudspeakers   = 65;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres


%Figure Output Settings
DocumentPath = 'tex\latex\APSIPA';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10);% + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 7/3;
FontSize = 9;
Font = 'Times';


%%

LUT_good = load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat'], ...
      'Quiet_SPL__Weight_Vs_Frequency', 'Frequencies');
  F_g = LUT_good.Frequencies;
LUT_good = LUT_good.Quiet_SPL__Weight_Vs_Frequency;

LUT_bad = load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '_unregulatedN.mat'], ...
      'Quiet_SPL__Weight_Vs_Frequency', 'Frequencies');
    F_b = LUT_bad.Frequencies;
LUT_bad = LUT_bad.Quiet_SPL__Weight_Vs_Frequency;


%%
% figure(1);
% plot(LUT_good(1,:),'b'); hold on;
% plot(LUT_good(end,:),'b');
% plot(LUT_bad(1,:),'r'); 
% plot(LUT_bad(end,:),'r'); hold off;

h=figure(1);
hold on;
set(gca, 'ColorOrderIndex',6); pl(3) = plot(F_g,min(LUT_good)-94,'-.', 'LineWidth', 0.5);
set(gca, 'ColorOrderIndex',6); pl(4) = plot(F_g,max(LUT_good)-94,'-','LineWidth', 0.5);
set(gca, 'ColorOrderIndex',7); pl(1) = plot(F_b,min(LUT_bad)-94,'-.', 'LineWidth', 0.5);
set(gca, 'ColorOrderIndex',7); pl(2) = plot(F_b,max(LUT_bad)-94,'-', 'LineWidth', 0.5);
 hold off;
%%
xlabel('Frequency (Hz)');
ylabel('SPL (dB)');
ylim([-70 0]);
grid on;box on;
set(gca,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1], ...
            'XScale', 'log');
        
%t(1).Position(2) = t(1).Position(2);
factor2 = 0.07;    %Page offset factor
factor1 = 2.3;  %Page height factor
offset_cm = FontSize*254/720 /10 * factor1;
set(gcf, 'PaperUnits','centimeters', ...
    'PaperSize', [plot_width, plot_width/aspect_ratio + offset_cm], ...
    'PaperPosition',[0 factor2 plot_width,plot_width/aspect_ratio + offset_cm]);

t(1) = title('Quiet Zone SPL Limits in Codebook', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
t(1).Position(2) = t(1).Position(2) + FontSize/2;

le = legend(pl([2 4]),{'{\it{N}} = 80','Regulated \it{N}'},'Location','southwest','Box','off', 'FontSize',FontSize, 'FontName',Font);
le.Position = le.Position + [0.043 -0.062 0 0];

print(['-d' print_fmt], [DocumentPath '\N_Planewave_LUT_Effect']);

% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;



