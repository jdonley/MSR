clc;
clear;
%close all;
tic;


%% Settings

% Video settings
writerObj = VideoWriter('example_Multitone.avi');
open(writerObj);
mov = [];

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
Resolution = 300;           % The higher the resolution the better the sampling rate and resultant image
NumberOfPlanewaves = [80 80];   % The more planewaves the more accurate the sound field

% Frequencies_in_zones = [ 500 1500 2000 2500 4000; ...
%                          300 1000 2200 3000 3500];
% Frequency_weights    = [ 0.2  1.5  2.5  5.0 10.0; ...
%                          0.1  1.2  2.8  7.0  8.0];  
                     
Frequencies_in_zones = [2000 ;...%4000; ...  % Zone 1 frequencies
                        1250 ];%3000];     % Zone 2 frequencies
Frequency_weights    = [   2.5   ;...%5 ; ...  % Importance weight for quiet zone of Zone 1
                           2.5   ];%5];      % Importance weight for quiet zone of Zone 2
                     
Angles_in_zones      = [15; ...         % Planewave angle Zone 1
                        90];           % Planewave angle Zone 2
Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 1.0; % Metres                    

[zones, frequencies] = size(Frequencies_in_zones);  % Currently only support for two zones using Orthogonal Basis Expansion

%zones = 1; % Force just one zone to be created. (Comment out for creation of both zones)


%%
repetitions = numden(sym(rat( max(Frequencies_in_zones) / min(Frequencies_in_zones) ))); %repetitions until phases line up again by finding least common multiple
repetitions=double(repetitions);
phase = pi/16:pi/16:2*pi*repetitions;

fprintf('\n====== Rendering ======\n\n');
fprintf('\tCompletion: ');n=0;
for p = 1:length(phase)
    sf = [];
    Soundfield = [];
    BrightField = [];
    
    for z = 1:zones
        for f = 1:frequencies
            %%
            phase_change = -phase(p)*Frequencies_in_zones(z, f)/max(Frequencies_in_zones);
            quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones(z, f), phase_change, Radius_of_zones(3-z), 'quiet');
            bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones(z, f), phase_change, Radius_of_zones(z),   'pw', 1.0, Angles_in_zones(z));
            quiet.res  = Resolution;
            bright.res = Resolution;
            quiet      = quiet.setDesiredSoundfield(true, 'suppress_output');
            if strcmp(DebugMode, 'DEBUG')
                bright = bright.setDesiredSoundfield(true, 'suppress_output');
            else
                bright = bright.setDesiredSoundfield(true, 'suppress_output');
            end
            
            %%
            sf = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
            sf = sf.addSpatialZone(quiet,  Position_of_zones(3-z,1), Position_of_zones(3-z,2));
            sf = sf.addSpatialZone(bright, Position_of_zones(z,1),   Position_of_zones(z,2));
            
            %%
            sf.QuietZ_Weight = Frequency_weights(z, f);
            sf = sf.setN(NumberOfPlanewaves(z));
            sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
            sf = sf.norm_soundfield;
            
            if (f==1) && (z==1)
                Soundfield = sf.Soundfield_desired;
                BrightField = sf.Bright_Field;
            else
                Soundfield = Soundfield + sf.Soundfield_desired;
                if (z==1)
                    BrightField = BrightField + sf.Bright_Field;
                end
            end
        end
    end
    tElapsed = toc;
    ratio = p / length(phase);
    tRem = (1-ratio) / ratio * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
    
    %%
    h=figure(1);
    sf.plotSoundfield(Soundfield);
    set(gcf,'Renderer','zbuffer');
    mov = [mov; getframe(h)];
    
end

%%
writeVideo(writerObj, repmat(mov,ceil(20/repetitions),1) );
close(writerObj);



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
