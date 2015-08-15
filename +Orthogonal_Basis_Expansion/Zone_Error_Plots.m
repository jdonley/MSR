clc;
clear;
close all;
tic;

%%
%Figure Output Settings
DocumentPath = 'tex\latex\APSIPA';
print_fmt = 'png'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10 + 6.35/10 + 88.9/10); %IEEE full page width
FontSize = 9;
Font = 'Times';
Adj= 254/720/10;

%%
Resolution = 400; % samples per metre

% Frequency = 2000; % Hz
% Weight = 1e4; % min of 1e-2, max of 1e4

Frequencies = [ 2000,  2000];   % Hz

Weights     = [  1e-2,  1e4]; % min of 1e-2, max of 1e4

%%
f_min = 150;
f_max = 8000;
N_min=28;
N_max = 300;

%%
h=figure(1);

%%
nam_val = {'FontSize',FontSize,'FontName',Font};
Rows = size(Weights,2);
Cols = size(Weights,1)*2; %Mag and Arg

for ind = 1:numel(Frequencies)
Frequency = Frequencies(ind);
Weight = Weights(ind);

N = ceil((N_max - N_min)*(f_max - f_min)^(-1.2)*(Frequency - f_min).^1.2 + N_min);

%%
%Soundfield settings
quiet  = Orthogonal_Basis_Expansion.spatial_zone(Frequency, 0, 0.30, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(Frequency, 0, 0.30, 'pw', 1.0, 15);
quiet.res  = Resolution;
bright.res = quiet.res;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');    
bright = bright.setDesiredSoundfield(true);

% phase = -rem( (0.3+0.3)*cos(bright.SourceOrigin.Angle), 343 / 2000) / (343 / 2000) * 360;
% temp = Orthogonal_Basis_Expansion.spatial_zone(2000, phase, 0.30, 'pw', 1.0, bright.SourceOrigin.Angle);
% temp.res = quiet.res;
% temp  =  temp.setDesiredSoundfield(true, 'suppress_output');   
% quiet.Soundfield_d = temp.Soundfield_d;
%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  0.60, 0);
soundfield = soundfield.addSpatialZone(bright, 0.60, 180);
%%
 soundfield.BrightZ_Weight     = 1.0;
 soundfield.QuietZ_Weight      = Weight;
 soundfield.UnattendedZ_Weight = 0.05;

soundfield = soundfield.setN(N);
soundfield = soundfield.createSoundfield('DEBUG', 1.0);

%%
%close all;
%h=figure(1);
%soundfield.plotSoundfield(soundfield.Soundfield_desired * exp( 1i * 0/180*pi),gray);

 

%saveas(h,'test.fig');
mask = ones(size(soundfield.Bright_Zone.Soundfield_d_mask));
mask(~soundfield.Bright_Zone.Soundfield_d_mask)=NaN;

actual = soundfield.Bright_Field;
desired = abs(soundfield.Bright_Zone.Soundfield_d).*exp(1i*(-angle(soundfield.Bright_Zone.Soundfield_d)));

plot_Size = Adj * 3.5;
y_off = Adj*3.0;
cb_x_off = -Adj * 0;

%
subplot(Rows,Cols,ind);
soundfield.plotSoundfield( (abs(desired) - abs(actual)).*mask, bone, false);

caxis([0 1]); grid minor;
cb=colorbar(nam_val{:}, 'Ticks',[0 0.25 0.5 0.75 1]);
set(gca,nam_val{:});
if (ind==1), t(ind)=title('\bf{Bright Zone $w_q=10^{-2}$}', nam_val{:},'interpreter','latex');
else t(ind)=title('\bf{Bright Zone $w_q=10^{4}$}', nam_val{:},'interpreter','latex'); end
xlabel('');%xlabel('Width (m)', nam_val{:});
if (ind==1),ylabel('Length (m)', nam_val{:});else ylabel(''); end
if (ind==1),colorbar off; end
pl_pos = get(gca,'Position'); 
if (ind==1), set(gca,'Position', pl_pos + [0 -y_off plot_Size plot_Size]);
else set(gca,'Position', pl_pos + [-0.07 -y_off plot_Size plot_Size]); end
if (ind~=1), pos=get(cb,'Position');set(cb,'Location','manual'); set(cb,'Position',pos + [cb_x_off 0.05 0 -0.1]);end
if (ind~=1), set(gca,'YTickLabel',''); end
set(gca,'XTickLabel','');

%
subplot(Rows,Cols,ind+2);
soundfield.plotSoundfield( (angle(desired ./ actual)).*mask, bone, false);

caxis([-pi pi]); grid minor;
cb=colorbar(nam_val{:}, 'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
set(gca,nam_val{:});title('');%t(ind+2)=title('Difference in Phase', nam_val{:});
xlabel('Width (m)', nam_val{:});
if (ind==1),ylabel('Length (m)', nam_val{:});else ylabel(''); end
if (ind==1),colorbar off; end
pl_pos = get(gca,'Position'); 
if (ind==1), set(gca,'Position', pl_pos + [0 0 plot_Size plot_Size]);
else set(gca,'Position', pl_pos + [-0.07 0 plot_Size plot_Size]); end
if (ind~=1), pos=get(cb,'Position');set(cb,'Location','manual'); set(cb,'Position',pos + [cb_x_off 0.05 0 -0.1]); end
if (ind~=1), set(gca,'YTickLabel',''); end
 
%%
fprintf('\tMSE in Bright zone: %.2fdB.\n', soundfield.Err_dB_Bright_Field);

fprintf('\tMSE in Quiet zone: %.2fdB.\n', soundfield.Err_dB_Quiet_Field);

fprintf('\tAcoustic Brightness Contrast: %.2fdB.\n', soundfield.Acoustic_Brightness_Contrast_dB);

fprintf('\n\tCalculated Angle in Bright Zone: %.1f°.\n', soundfield.Angle_Bright_Field);

fprintf('\tCalculated Angle in Quiet Zone: %.1f°.\n', soundfield.Angle_Quiet_Field);


end

%%

paper_size_x = 0.0;
paper_size_y = -10.0;

fig_size_x = 0.0;
fig_size_y = 0.0;

offset_x = 0.0;
offset_y = 0.0;


paperAdj_cm_x =  Adj * paper_size_x;
paperAdj_cm_y = Adj * paper_size_y;
figAdj_cm_x = Adj * fig_size_x;
figAdj_cm_y = Adj * fig_size_y;
offset_cm_x = Adj * offset_x;
offset_cm_y = Adj * offset_y;

title_offset = FontSize*Adj *10;

for i=1:numel(Frequencies), t(i).Position(2) = t(i).Position(2)+title_offset;end
t(2).Position(1) = t(2).Position(1)+6.5;
%tightfig(h);

set(gcf, 'PaperUnits','centimeters', ...
         'PaperSize',       [plot_width/Rows + paperAdj_cm_x, ...
                            plot_width/Cols + paperAdj_cm_y], ...
         'PaperPosition',   [offset_cm_x, offset_cm_y, ...
                            plot_width/Rows + figAdj_cm_x, ...
                            plot_width/Cols + figAdj_cm_y]);
 set(gcf, 'Units','centimeters', ...
         'OuterPosition',   [offset_cm_x, offset_cm_y, ...
                            plot_width/Rows + figAdj_cm_x, ...
                            plot_width/Cols + figAdj_cm_y]);                 

print(h, ['-d' print_fmt], ['-r' num2str(print_res)], [DocumentPath '\Spatial_Error_Improvement']);
%print(h, ['-d' print_fmt], [DocumentPath '\Spatial_Error_Improvement']);
close(h);
% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script


%%
