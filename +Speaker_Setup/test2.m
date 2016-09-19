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

array_type = '2line';
spkr_radius = 1.3;
spkr=1;

if strcmp(array_type, 'circle')
    [x,y] = pol2cart(-90/180*pi, 0.6);
    x_ = sqrt(spkr_radius^2-y^2);
    th_c = atan2(y,-x_);
    th = th_c;
    spkr_spacing = []; %Auto-calculate spacing
elseif ~isempty(strfind(array_type, 'line'))
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
f = 1000;
c=343;
setup = Speaker_Setup.createSetup({...
    'frequency',                    f, ...
    speech_layout{:}, ...
    array{:}, ...
    'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
    'reproduction_radius',          1.0, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 1e4, ...
    'unattended_weight',            0.05, ...
    'brightzone_radius',            0.30, ...
    'brightzone_source_dist',       x_, ...
    'brightzone_pos_distance',      0.6, ...
    'quietzone_radius',             0.3, ...
    'quietzone_pos_distance',       0.6, ...
    'maximum_frequency',            8000, ...
    'angleof_loudspeakerarrcentre', 180, ...
    'loudspeaker_object',           Parametric_Synthesis.parametric_soundfield});

setup.ExtendedField = true;
setup.Origin   = [0.0,1.3];
%%
setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');

setup = setup.calc_Loudspeaker_Weights();
[x,y]=pol2cart(setup.Loudspeaker_Locations([1 end],1),setup.Loudspeaker_Locations([1 end],2));
d = sum(diff([x,y]).^2).^0.5;
setup.Loudspeaker_Weights(1:end/2) = setup.Loudspeaker_Weights(1:end/2)*exp(1i* (-d*f/c * pi+pi));
setup = setup.reproduceSoundfield('DEBUG');

%%
figNums = [101,102];
hold off;
realistic = false;
details.DrawDetails = false;
details.zoneLineWid = 1.5;
details.arrowLineWid = 0.4;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 30;

%pk = max(abs((setup.Bright_Samples(:))));
pk = max(abs((setup.Quiet_Samples(:))));
Z = setup.Soundfield_reproduced;
Z2 = abs(Z);

Z_ = mag2db((Z)./pk);
figure(figNums(1))
setup.plotSoundfield( Z, 'scientific_D1A', realistic, details);
figure(figNums(2))
setup.plotSoundfield( Z2, 'scientific_L9', realistic, details);

% for fn = 1:numel(figNums)
%     figure(figNums(fn));
%     if spkr ==1
%         R = [0 2].*spkr_radius*100 ; xlim(R);ylim(R);
%     elseif spkr == 2
%         R = [0 2].*spkr_radius*100 - (spkr_radius-x_)*100; xlim(R);ylim(R);
%     end
% end
%  caxis([-30, 0] );


%title('Small Zone Weight');

%  hold on;

%  Z = setup.Soundfield_reproduced;
%  %Z_ = -QualityGuidedUnwrap2D_r1(Z);
%  Z_ = GoldsteinUnwrap2D_r1(Z);
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
%  quiver3( X , Y, ones(size(U))*4, U .* abs(Zs(2:end-1,2:end-1)) , V .* abs(Zs(2:end-1,2:end-1)), zeros(size(U)), 1, 'k' );
%quiver3( X , Y, ones(size(U))*4, U  , V , zeros(size(U)), 1, 'k' );


disp(['Contrast: ' num2str(mag2db(setup.Acoustic_Contrast)) 'dB']);
disp(['     MSE: ' num2str(mag2db(setup.MSE_Bright)) 'dB']);

%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script