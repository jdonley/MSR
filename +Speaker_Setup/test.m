clc;
%clear;
%close all;
tic;

%%
f_max = 8000;
quiet  = Orthogonal_Basis_Expansion.spatial_zone(6000, 0, 0.3, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(6000, 0, 0.3, 'pw', 1.0, 20);
quiet.res  = 50;
bright.res = 50;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');    
bright = bright.setDesiredSoundfield(true);

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  0.6, 0);
soundfield = soundfield.addSpatialZone(bright, 0.6, 180);
%%
 soundfield.BrightZ_Weight     = 1.0;
 soundfield.QuietZ_Weight      = 10;
 soundfield.UnattendedZ_Weight = 0.05;

soundfield = soundfield.setN(230);
Radius = 1.0;
soundfield = soundfield.createSoundfield('DEBUG', Radius);

%%
%figure(1);
%soundfield.plotSoundfield();


%%
 setup = Speaker_Setup.loudspeaker_setup;
 setup = setup.addMultizone_Soundfield(soundfield);
 %setup.Loudspeaker_Count = 2*ceil(f_max/343 * exp(1) * Radius / 2 ) + 1 ;%16;
 setup.Speaker_Arc_Angle = 180;
 setup.Angle_FirstSpeaker = 120;
 k = f_max/343*2*pi;
 M = ceil(k*Radius);
 setup.Loudspeaker_Count = 32;%ceil( setup.Speaker_Arc_Angle/360 * 2*M + 1 ) ;%16;
 setup = setup.setRadius(1.5);
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');
 %%
 figure(1);
 setup.plotSoundfield(setup.Soundfield_reproduced);

 %  figure(3);
%  plot(abs(setup.Loudspeaker_Weights));
%  figure(4);
%  plot(angle(setup.Loudspeaker_Weights)/pi*180);
 
 
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script