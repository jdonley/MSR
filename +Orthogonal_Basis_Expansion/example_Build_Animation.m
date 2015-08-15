clc;
%clear;
%close all;
tic;

%%
writerObj = VideoWriter('example_Monotone.avi');
open(writerObj);
phase = pi/16:pi/16:2*pi;
mov = [];


%%
for p = 1:length(phase)
    
    quiet  = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'quiet');
    quiet.res  = 200;
    quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');
    bright = Orthogonal_Basis_Expansion.spatial_zone(2000, -phase(p), 0.30, 'pw', 1.0, 15);
    bright.res = 200;
    bright = bright.setDesiredSoundfield(true);
    
    %%
    soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
    soundfield = soundfield.addSpatialZone(quiet,  0.60, 0);
    soundfield = soundfield.addSpatialZone(bright, 0.60, 180);
    
    %%
    soundfield = soundfield.createSoundfield('DEBUG', 1.0);
    
    %%
    h=figure(1);
    soundfield.plotSoundfield();
    set(gcf,'Renderer','zbuffer');
    mov = [mov; getframe(h)];
    
end

%%
writeVideo(writerObj, repmat(mov,20,1) );
close(writerObj);



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
