clc;
clear;
%close all;
tic;

%%

        
f = 3500;
w = 1e4;
col = 'k';
f_max = 8000;

quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.30, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.30, 'pw', 1.0, 15);
quiet.res  = 100;
bright.res = quiet.res;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');    
bright = bright.setDesiredSoundfield(true);

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  0.6, 0);
soundfield = soundfield.addSpatialZone(bright, 0.6, 180);
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
 
 k = f_max/343*2*pi;
 M = ceil(k*Radius);
setup.Loudspeaker_Count = 24; %ceil( setup.Speaker_Arc_Angle/360 * 2*M + 1 ) %16;
setup = setup.setRadius(1.5);

setup.Speaker_Arc_Angle = 180;
setup.Angle_FirstSpeaker = 90;

setup = setup.setLoudspeakerType('Genelec 8010A');
setup.Speaker_Array_Centre = 180;
setup = setup.setLoudspeakerSpacing( 0.01 );
 
 spkr_first = [-1.49, -1.49];
 spkr_last = [-1.49, 1.49];
 x = linspace(spkr_first(1), spkr_last(1),setup.Loudspeaker_Count);
 y = linspace(spkr_first(2), spkr_last(2),setup.Loudspeaker_Count);
 y = ((0.121+0.005)/2):(0.121+0.005):(0.121+0.005)*setup.Loudspeaker_Count/2;
 y = [flip(-y), y];
 [spkr_theta, spkr_rho]=cart2pol(x,y);
 setup.Loudspeaker_Locations = [spkr_theta', spkr_rho'];
 
 %setup = setup.calc_Loudspeaker_Locations();
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

%  Field{i} = abs(setup.Soundfield_reproduced);
%  magdif{i,j} = setup.Bright_Sample - setup.Quiet_Sample;
% end
%  end
%  end
% F = zeros(size(Field{1}));
% for i = 1:length(Field)
%    F = F+Field{i};
% end
% F = F/length(Field);
 %%
 figure(1);
 realistic = false;
 setup.plotSoundfield( setup.Soundfield_reproduced, 'default', realistic);
 %caxis([0 1.1]);
 %mean([magdif{:}])
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
export_fig('test.png',gcf);
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script