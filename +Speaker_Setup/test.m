clc;
clear;
%close all;
tic;

%%
i=0;
j=0;
 for s = 1:1
 ang=[-30, 30];
 for f=logspace(log10(150),log10(400),5)
i=i+1;
for w=logspace(log10(1),log10(1e3),5)
j=j+1;
%f = 200;
col = 'k';
f_max = 8000;

quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.30, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.30, 'pw', 1.0, ang(s));
quiet.res  = 50;
bright.res = quiet.res;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');    
bright = bright.setDesiredSoundfield(true);

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  0.6, 180);
soundfield = soundfield.addSpatialZone(bright, 0.6, 0);
%%
 soundfield.BrightZ_Weight     = 1.0;
 soundfield.QuietZ_Weight      = w;
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
 setup.Speaker_Arc_Angle = 180;
 setup.Angle_FirstSpeaker = 90;
 k = f_max/343*2*pi;
 M = ceil(k*Radius);
 setup.Loudspeaker_Count = 32; %ceil( setup.Speaker_Arc_Angle/360 * 2*M + 1 ) %16;
 setup = setup.setRadius(1.5);
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

 Field{i} = abs(setup.Soundfield_reproduced);
 magdif{i,j} = setup.Bright_Sample - setup.Quiet_Sample;
end
 end
 end
F = zeros(size(Field{1}));
for i = 1:length(Field)
   F = F+Field{i};
end
F = F/length(Field);
 %%
 figure(1);
 realistic = false;
 setup.plotSoundfield( F, 'default', realistic);
 caxis([0 1.1]);
 mean([magdif{:}])
% figure(2);
% plot(abs(setup.Loudspeaker_Weights), col);hold on;
% figure(3);
% scatter(f,mean(abs(setup.Bright_Samples(:))), col); hold on;
% scatter(f,mean(abs(setup.Bright_Samples(:)))+ std(abs(setup.Bright_Samples(:))), 'r'); 
% scatter(f,mean(abs(setup.Bright_Samples(:)))- std(abs(setup.Bright_Samples(:))), 'r'); 
% drawnow();
%  figure(3);
%  plot(abs(setup.Loudspeaker_Weights));
%  figure(4);
%  plot(angle(setup.Loudspeaker_Weights)/pi*180);
 
%end 
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script