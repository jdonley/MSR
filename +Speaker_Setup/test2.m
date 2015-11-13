clc;
clear;
%close all;
tic;

%%
room_width = 3.200; %m
rad = ( room_width - 2*0.115 - 2*0.085 ) / 2;


% speech_layout = {'brightzone_pos_angle',        180, ...
%                  'quietzone_pos_angle',         0, ...
%                  'brightzone_source_angle',     15};
% masker_layout = {'brightzone_pos_angle',        0, ...
%                  'quietzone_pos_angle',         180, ...
%                  'brightzone_source_angle',     0};

speech_layout = {'brightzone_pos_angle',        90, ...
                 'quietzone_pos_angle',         -90, ...
                 'brightzone_source_angle',     0};
masker_layout = {'brightzone_pos_angle',        -90, ...
                 'quietzone_pos_angle',         90, ...
                 'brightzone_source_angle',     0};
             
setup = Speaker_Setup.createSetup({...
            'frequency',                    2000, ...
            speech_layout{:}, ...
            'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
            'reproduction_radius',          1.0, ...
            'bright_weight',                1.0, ...
            'quiet_weight',                 1e4, ...
            'unattended_weight',            0.05, ...
            'brightzone_radius',            0.3, ...
            'brightzone_pos_distance',      0.6, ...
            'quietzone_radius',             0.3, ...
            'quietzone_pos_distance',       0.6, ...
            'numberof_loudspeakers',        24, ...
            'loudspeaker_radius',           1.5, ...
            'maximum_frequency',            8000, ...
            'angleto_firstloudspeaker',     90, ...
            'angleof_loudspeakerarc',       180 + 180/24, ...
            'loudspeaker_model',            'Genelec 8010A', ...
            'angleof_loudspeakerarrcentre', 180, ...
            'loudspeaker_spacing',          0.01, ...
            'speaker_array_type',           'line'});
        
%%
setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

 
 %%
 hold off;
 figure(1);
 realistic = false;
 setup.plotSoundfield( abs(setup.Soundfield_reproduced), 'default', realistic);
 
 hold on;
 
 Z = setup.Soundfield_reproduced;
 %Z_ = -QualityGuidedUnwrap2D_r1(Z);
 Z_ = -GoldsteinUnwrap2D_r1(Z);
 
 Zs = abs(Z);
 Zs = Zs(1:10:end,1:10:end);
 Z__ = Z_(1:10:end,1:10:end);
 
 [U,V] = gradient( Z__(2:end-1,2:end-1) );
 [X,Y] = meshgrid( 1:size(Z,1) , 1:size(Z,2) );
 X = X(1:10:end,1:10:end);
 Y = Y(1:10:end,1:10:end);
 X = X(2:end-1,2:end-1);
 Y = Y(2:end-1,2:end-1);
 quiver3( X , Y, ones(size(U))*10, U .* abs(Zs(2:end-1,2:end-1)) , V .* abs(Zs(2:end-1,2:end-1)), zeros(size(U)), 1, 'k' );
 %quiver3( X , Y, ones(size(U))*10, U  , V , zeros(size(U)), 1, 'k' );
 
%  figure(1);
%  for p = 0:10:1080
% Z = (setup.Soundfield_reproduced .* exp(1i*p/180*pi));
% hold off;
% setup.plotSoundfield( abs(Z)/max(abs(Z(:))) .* angle(Z), 'default', realistic); hold on;
%  quiver3( X , Y, ones(size(U))*10, U .* abs(Zs(2:end-1,2:end-1)) , V .* abs(Zs(2:end-1,2:end-1)), zeros(size(U)), 1, 'k' );
% drawnow();
% end

disp(['Contrast: ' num2str((setup.Bright_Sample - setup.Quiet_Sample)*100) '%']);

%%
fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script