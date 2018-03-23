clc;
clear;
close all;
% N Planewave effect on LUT Comparison

LUT_resolution = '512f_256w';
angle_pw = 15;
loudspeakers   = 295;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres


%Figure Output Settings
DocumentPath = 'tex\latex\ChinaSIP';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10);% + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 7/3;
FontSize = 9;
Font = 'Times';

%%
LUT = load(['Z:\+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat'], ...
      'Quiet_SPL__Weight_Vs_Frequency', 'Frequencies', 'Weights');
  F = LUT.Frequencies;
  W = LUT.Weights;
LUT = LUT.Quiet_SPL__Weight_Vs_Frequency - 94;


%%
surf(F,W,LUT,'LineStyle','none');h = gcf;view(2);
colormap bone;
c = colorbar;
ax=gca;
set(ax,'XScale','log','YScale','log');

t = title('SPL in Quiet Zone');

xlabel('Frequency (Hz)');
ylabel('Weight');
ylabel(c, 'SPL (dB)');

ylim([min(W) max(W)]);
%xlim([min(F) max(F)]);
c.Ticks = [-55 -30 -5];
ax.YTick = [1e-2 1e1 1e4];

grid on;box on;
set(ax,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1]);
        
% Colorbar width
% axpos = ax.Position;
% cpos = c.Position;
% cpos(3) = cpos(3) * 0.5;
% c.Position = cpos;
% ax.Position = axpos;

%%
factor1 = -0.4;  %Page height factor
factor2 = 0.8;  %Page height factor
factor3 = 0.0;    %Page offset factor
offset1_cm = FontSize*254/720 /10 * factor1;
offset2_cm = FontSize*254/720 /10 * factor2;
set(h, 'PaperUnits','centimeters', ...
    'PaperSize', [plot_width, plot_width/aspect_ratio + offset1_cm], ...
    'PaperPosition',[0 factor3 plot_width, plot_width/aspect_ratio + offset2_cm]);

t.Position(2) = t.Position(2)*2;

%print(['-d' print_fmt], ['-r' num2str(print_res)], [DocumentPath '\Example_LUT']);
% print(['-d' print_fmt], ['-r' num2str(print_res)], [DocumentPath '\Example_High_Res_LUT']);

close(h);

% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;