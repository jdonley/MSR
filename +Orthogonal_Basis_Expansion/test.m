%clc;
clear;
%close all;
tic;

%%
%Figure Output Settings
DocumentPath = 'tex\latex\APSIPA';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10 + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 5/3;
FontSize = 9;
Font = 'Times';

%%
%Soundfield settings
quiet  = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'pw', 1.0, 0);
quiet.res  = 50;
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
 soundfield.QuietZ_Weight      = 2.5;
 soundfield.UnattendedZ_Weight = 0.05;

soundfield = soundfield.setN(109);
soundfield = soundfield.createSoundfield('DEBUG', 1.0);

%%
%close all;
h=figure(1);
soundfield.plotSoundfield(abs(soundfield.Soundfield_desired * exp( 1i * 0/180*pi)),parula);

 

%saveas(h,'test.fig');
% figure(3)
% soundfield.plotSoundfield( soundfield.Quiet_Error_Field );

fprintf('\tMSE in Bright zone: %.2fdB.\n', soundfield.Err_dB_Bright_Field);

fprintf('\tMSE in Quiet zone: %.2fdB.\n', soundfield.Err_dB_Quiet_Field);

fprintf('\tAcoustic Brightness Contrast: %.2fdB.\n', soundfield.Acoustic_Brightness_Contrast_dB);

fprintf('\n\tCalculated Angle in Bright Zone: %.1f°.\n', soundfield.Angle_Bright_Field);

fprintf('\tCalculated Angle in Quiet Zone: %.1f°.\n', soundfield.Angle_Quiet_Field);


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script


%%
