clc;
%clear;
%close all;
tic;

%%
for f=logspace(log10(150),log10(2000),40)
col = 'k';
f_max = 8000;
%f = 450;
quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.3, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.3, 'pw', 1.0, 0);
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
 soundfield.QuietZ_Weight      = 0.01;
 soundfield.UnattendedZ_Weight = 0.05;

[N, Frequencies] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( 512, 28, 300, 150, f_max );
[~,I]=find(Frequencies >= f);
soundfield = soundfield.setN( N(I(1)) );
Radius = 1.0;
soundfield = soundfield.createSoundfield('DEBUG', Radius);

%%
%figure(1);
%soundfield.plotSoundfield();


%%
 setup = Speaker_Setup.loudspeaker_setup;
 setup = setup.addMultizone_Soundfield(soundfield);
 %setup.Loudspeaker_Count = 2*ceil(f_max/343 * exp(1) * Radius / 2 ) + 1 ;%16;
 setup.Speaker_Arc_Angle = 360;
 setup.Angle_FirstSpeaker = 0;
 k = f_max/343*2*pi;
 M = ceil(k*Radius);
 setup.Loudspeaker_Count = ceil( setup.Speaker_Arc_Angle/360 * 2*M + 1 ) ;%16;
 setup = setup.setRadius(1.5);
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

 %%
 figure(1);
 setup.plotSoundfield(abs(setup.Soundfield_reproduced));
figure(2);
plot(abs(setup.Loudspeaker_Weights), col);hold on;
figure(3);
scatter(f,mean(abs(setup.Bright_Samples(:))), col); hold on;
scatter(f,mean(abs(setup.Bright_Samples(:)))+ std(abs(setup.Bright_Samples(:))), 'r'); 
scatter(f,mean(abs(setup.Bright_Samples(:)))- std(abs(setup.Bright_Samples(:))), 'r'); 
drawnow();
%  figure(3);
%  plot(abs(setup.Loudspeaker_Weights));
%  figure(4);
%  plot(angle(setup.Loudspeaker_Weights)/pi*180);
 
end 
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script