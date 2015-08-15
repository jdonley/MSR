clc;
%clear;
close all;
tic;

%%
writerObj = VideoWriter('example_Varied_Weighting.avi');
open(writerObj);
resolution = 300;
f=2000;
QZ_pos_a = 0;
weights = [logspace(log10(1e-2),log10(1e4),200) logspace(log10(1e4),log10(1e-2),200)];
mov = [];

figure('units','points','position',[0 0 1000 1000]);

%%
fprintf('\n====== Rendering ======\n\n');
fprintf('\tCompletion: ');n=0;
p_last=0;
for index = 1:length(weights)
    
    W = weights(index);
    
    quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.30, 'quiet');
    quiet.res  = resolution;
    quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');
    bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, 0.30, 'pw', 1.0, 15);
    bright.res = resolution;
    bright = bright.setDesiredSoundfield(true, 'suppress_output');
    
    %%
    soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
    soundfield = soundfield.addSpatialZone(quiet,  0.60, QZ_pos_a);
    soundfield = soundfield.addSpatialZone(bright, 0.60, 180);
    
    %% 
    soundfield.BrightZ_Weight     = 1.0;
    soundfield.QuietZ_Weight      = W;
    soundfield.UnattendedZ_Weight = 0.05;
    soundfield = soundfield.setN(80);
    soundfield = soundfield.createSoundfield('DEBUG', 1.0);
    
    %%
    setup = Speaker_Setup.loudspeaker_setup;
    setup = setup.addMultizone_Soundfield(soundfield);
    setup.Loudspeaker_Count = 32;
    setup.Speaker_Arc_Angle = 180;
    setup.Angle_FirstSpeaker = 90;
    setup = setup.setRadius(1.5);
    
    setup = setup.calc_Loudspeaker_Weights();
    setup = setup.reproduceSoundfield('DEBUG');
    
    %%    
    for p = (p_last+20):20:(p_last+40)
        h=figure(1);
        setup.plotSoundfield( setup.Soundfield_reproduced * exp(-1i*p/180*pi) );
        title('Multizone Soundfield Occlusion Problem');
        set(gcf,'Renderer','zbuffer');
        mov = [mov; getframe(h)];
    end
    p_last = p;

    %%
    tElapsed = toc;
    ratio = p / length(weights);
    tRem = (1-ratio) / ratio * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
    
end

%%
writeVideo(writerObj, mov );
close(writerObj);



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
