%clc;
clear;
%close all;
tic;

%%
room_width = 3.200; % metres
rad = ( room_width - 2*0.115 - 2*0.085 ) / 2; % metres

% speech_layout = {'brightzone_pos_angle',        180, ...
%                  'quietzone_pos_angle',         0, ...
%                  'brightzone_source_angle',     14.5};
% masker_layout = {'brightzone_pos_angle',        0, ...
%                  'quietzone_pos_angle',         180, ...
%                  'brightzone_source_angle',     0};

% speech_layout = {'brightzone_pos_angle',        90, ...
%                  'quietzone_pos_angle',         -90, ...
%                  'brightzone_source_angle',     0};
% masker_layout = {'brightzone_pos_angle',        -90, ...
%                  'quietzone_pos_angle',         90, ...
%                  'brightzone_source_angle',     0};
             
speech_layout = {'brightzone_pos_angle',        90, ...
                 'quietzone_pos_angle',         -90, ...
                 'brightzone_source_angle',     0, ...    
                 'brightzone_source_type',      'pw'};
             
masker_layout = {'brightzone_pos_angle',        -90, ...
                 'quietzone_pos_angle',         90, ...
                 'brightzone_source_angle',     180, ...    
                 'brightzone_source_type',      'ps'};

array_type = 'line';
spkr_radius = 1.3;
spkr=1;

if strcmp(array_type, 'circle')
    [x,y] = pol2cart(-90/180*pi, 0.6);
    x_ = sqrt(spkr_radius^2-y^2);
    th_c = atan2(y,-x_);
    th = th_c;
    spkr_spacing = []; %Auto-calculate spacing
elseif strcmp(array_type, 'line')
    x_=spkr_radius;
    th_c = 180;
    th = atan2(-0.6,-spkr_radius);    
    spkr_spacing = 0.001; %1mm spacing between adjacent loudspeakers
end
if spkr ==1
    %regular
    N_spkrs = 24;
    array = { ...
        'angleto_firstloudspeaker',      90, ...
        'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
        'numberof_loudspeakers',         N_spkrs, ...
        'loudspeaker_model',             'Genelec 8010A', ...
        'loudspeaker_radius',            spkr_radius, ...
        'loudspeaker_spacing',           spkr_spacing, ...
        'speaker_array_type',            array_type};
elseif spkr == 2
    %parametric
    array = { ...
        'angleto_firstloudspeaker',     th/pi*180, ...
        'loudspeaker_radius',           x_, ...
        'numberof_loudspeakers',        1, ...
        'loudspeaker_model',            'Parametric', ...
        'loudspeaker_spacing',          0.01, ...
        'speaker_array_type',           'line'};
end
             
setup = Speaker_Setup.createSetup({...
            'frequency',                    1000, ...
            masker_layout{:}, ...
            array{:}, ...
            'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
            'reproduction_radius',          1.0, ...
            'bright_weight',                1.0, ...
            'quiet_weight',                 0, ...
            'unattended_weight',            0.05, ...
            'brightzone_radius',            0.3, ...
            'brightzone_source_dist',       x_, ...
            'brightzone_pos_distance',      0.6, ...
            'quietzone_radius',             0.3, ...
            'quietzone_pos_distance',       0.6, ...
            'maximum_frequency',            8000, ...
            'angleof_loudspeakerarrcentre', 180, ...
            'loudspeaker_object',           Parametric_Synthesis.parametric_soundfield});
        
%%
setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

 %%
 figure(123);
 hold off;
 realistic = false;
 
 pk = max(abs(real(setup.Bright_Samples(:))));
 Z = setup.Soundfield_reproduced;
 
 setup.plotSoundfield( (Z), 'default', realistic);
 
if spkr ==1
    R = [0 2].*spkr_radius*100 ; xlim(R);ylim(R);
elseif spkr == 2
    R = [0 2].*spkr_radius*100 - (spkr_radius-x_)*100; xlim(R);ylim(R);
end
 %caxis([-pk, pk] );
 
 
 %title('Small Zone Weight');
 
 hold on;
 
 Z = setup.Soundfield_reproduced;
 %Z_ = -QualityGuidedUnwrap2D_r1(Z);
 Z_ = GoldsteinUnwrap2D_r1(Z);
 
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


disp(['Contrast: ' num2str((setup.Bright_Sample - setup.Quiet_Sample)*100) '%']);
disp(['          ' num2str(mag2db(setup.Acoustic_Contrast)) 'dB']);

%%
fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script