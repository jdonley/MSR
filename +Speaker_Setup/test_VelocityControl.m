clc;
%clear;
close all;
tic;

%%
f_max = 8000;
quiet  = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'pw', 1.0, 15 + 180);
quiet.res  = 100;
bright.res = 100;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');    
bright = bright.setDesiredSoundfield(true);

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  0.60, 0);
soundfield = soundfield.addSpatialZone(bright, 0.60, 180);
%%
 soundfield.BrightZ_Weight     = 1.0;
 soundfield.QuietZ_Weight      = 2.5;
 soundfield.UnattendedZ_Weight = 0.05;

soundfield = soundfield.setN(80);
Radius = 1.0;
soundfield = soundfield.createSoundfield('DEBUG', Radius);

%%
%figure(1);
%soundfield.plotSoundfield();


%%
 setup = Speaker_Setup.loudspeaker_setup_VelocityControl;
 setup = setup.addMultizone_Soundfield(soundfield);
 setup.Loudspeaker_Count = 2*ceil(f_max/343 * exp(1) * Radius / 2) + 1;%16;
 setup.Speaker_Arc_Angle = 360;
 setup.Angle_FirstSpeaker = 105;
 setup = setup.setRadius(1.5);
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');
 %%
 hold off;
 p  = setup.Soundfield_PressureReproduced; %Pressure
 v_ = setup.Soundfield_VelocityReproduced;  %Velocity
 
 figure(1);
 setup.plotSoundfield(real(p));
 
 figure(2);
 setup.plotSoundfield(abs(v_));
 
 figure(3);
 I = 1/2 * real( p .* conj(v_) );
 setup.plotSoundfield( I );
 
% figure
% [DX, DY] = gradient(-abs(setup.Soundfield_VelocityReproduced),1/setup.res,1/setup.res);
% w = size(setup.Soundfield_VelocityReproduced,1);
% [X,Y] = meshgrid(linspace(-w/2/setup.res,w/2/setup.res,w));
% [sx,sy] = meshgrid(linspace(-w/2/setup.res,w/2/setup.res,w*0.5));
% streamline(stream2(X,Y,DX,DY,sx,sy,[0.1 10000]));
%  figure(3);
%  plot(abs(setup.Loudspeaker_Weights));
%  figure(4);
%  plot(angle(setup.Loudspeaker_Weights)/pi*180);
 
 
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script