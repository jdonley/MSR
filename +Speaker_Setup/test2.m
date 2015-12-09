%clc;
clear;
%close all;
tic;

%%
room_width = 3.200; % metres
rad = ( room_width - 2*0.115 - 2*0.085 ) / 2; % metres


speech_layout = {'brightzone_pos_angle',        180, ...
                 'quietzone_pos_angle',         0, ...
                 'brightzone_source_angle',     14.5};
masker_layout = {'brightzone_pos_angle',        0, ...
                 'quietzone_pos_angle',         180, ...
                 'brightzone_source_angle',     0};

% speech_layout = {'brightzone_pos_angle',        90, ...
%                  'quietzone_pos_angle',         -90, ...
%                  'brightzone_source_angle',     0};
% masker_layout = {'brightzone_pos_angle',        -90, ...
%                  'quietzone_pos_angle',         90, ...
%                  'brightzone_source_angle',     -60};
             
setup = Speaker_Setup.createSetup({...
            'frequency',                    2000, ...
            speech_layout{:}, ...
            'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
            'reproduction_radius',          1.0, ...
            'bright_weight',                1.0, ...
            'quiet_weight',                 0, ...
            'unattended_weight',            0.05, ...
            'brightzone_radius',            0.3, ...
            'brightzone_pos_distance',      0.6, ...
            'quietzone_radius',             0.3, ...
            'quietzone_pos_distance',       0.6, ...
            'numberof_loudspeakers',        65, ...
            'loudspeaker_radius',           1.5, ...
            'maximum_frequency',            8000, ...
            'angleto_firstloudspeaker',     90, ...
            'angleof_loudspeakerarc',       360 , ...
            'loudspeaker_model',            'Genelec 8010A', ...
            'angleof_loudspeakerarrcentre', 180, ...
            'loudspeaker_spacing',          [], ...
            'speaker_array_type',           'circle'});
        
%%
setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

 %%
 hold off;
 figure(123);
 realistic = false;
 
 pk = max(abs(real(setup.Bright_Samples(:))));
 Z = setup.Soundfield_reproduced;
 
 setup.plotSoundfield( (Z), 'default', realistic);
 
 caxis([-pk, pk] );
 title('Small Zone Weight');
 
%  hold on;
%  
%  Z = setup.Soundfield_reproduced;
%  %Z_ = -QualityGuidedUnwrap2D_r1(Z);
%  Z_ = -GoldsteinUnwrap2D_r1(Z);
%  
%  Zs = abs(Z);
%  Zs = Zs(1:10:end,1:10:end);
%  Z__ = Z_(1:10:end,1:10:end);
%  
%  [U,V] = gradient( Z__(2:end-1,2:end-1) );
%  [X,Y] = meshgrid( 1:size(Z,1) , 1:size(Z,2) );
%  X = X(1:10:end,1:10:end);
%  Y = Y(1:10:end,1:10:end);
%  X = X(2:end-1,2:end-1);
%  Y = Y(2:end-1,2:end-1);
%  quiver3( X , Y, ones(size(U))*10, U .* abs(Zs(2:end-1,2:end-1)) , V .* abs(Zs(2:end-1,2:end-1)), zeros(size(U)), 1, 'k' );
%  %quiver3( X , Y, ones(size(U))*10, U  , V , zeros(size(U)), 1, 'k' );


disp(['Contrast: ' num2str((setup.Bright_Sample - setup.Quiet_Sample)*100) '%']);
disp(['          ' num2str(mag2db(setup.Acoustic_Contrast)) 'dB']);

%%
fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script