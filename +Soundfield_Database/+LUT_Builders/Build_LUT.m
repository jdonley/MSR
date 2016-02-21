clc;
%clear;
close all;
tic; %Start timing this script
delete(gcp);
current_pool = parpool; %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))



%% Settings
% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
% You can suppress the output when executing a complete soundfield construction with 'suppress_output'
%%%% OLD LAYOUT %%%%
% speech_layout = {'brightzone_pos_angle',        180, ...
%                  'quietzone_pos_angle',         0, ...
%                  'brightzone_source_angle',     15};
% masker_layout = {'brightzone_pos_angle',        0, ...
%                  'quietzone_pos_angle',         180, ...
%                  'brightzone_source_angle',     0};

%%%% NEW LAYOUT %%%%
%LUTlayouts = [1, 2, 3];
LUTlayouts = [3];
for LUTBUILD = LUTlayouts
    
    % Multizone Speech
    if LUTBUILD == 1
        spkr = 1;
        zonelayout = 'speech';
    % Multizone Masker
    elseif LUTBUILD == 2
        spkr = 1;
        zonelayout = 'masker';
    % Parametric Masker
    elseif LUTBUILD == 3
        spkr = 2;
        zonelayout = 'masker';
    end

array_type = 'circle';
    
if strcmp(zonelayout,'speech')
    layout = {'brightzone_pos_angle',        90, ...
              'quietzone_pos_angle',         -90, ...
              'brightzone_source_angle',     0, ...    
              'brightzone_source_type',      'pw'};

elseif strcmp(zonelayout,'masker')
    layout = {'brightzone_pos_angle',        -90, ...
              'quietzone_pos_angle',         90, ...
              'brightzone_source_angle',     180, ...
              'brightzone_source_type',      'ps'};
end

if strcmp(array_type, 'circle')
    [x,y] = pol2cart(-90/180*pi, 0.6);
    x_ = sqrt(1.5^2-y^2);
    th_c = atan2(y,-x_);
    th = th_c;
elseif strcmp(array_type, 'line')
    x_=1.5;
    th_c = 180;
    th = atan2(-0.6,-1.5);
end
if spkr ==1
    %regular
    array = {'angleto_firstloudspeaker',      90, ...
        'numberof_loudspeakers',         24, ...
        'loudspeaker_model',             'Genelec 8010A', ...
        'loudspeaker_radius',            1.5, ...
        'loudspeaker_spacing',           [], ...
        'speaker_array_type',            array_type};
elseif spkr == 2
    %parametric
    array = {'angleto_firstloudspeaker',     th/pi*180, ...
        'loudspeaker_radius',           x_, ...
        'numberof_loudspeakers',        1, ...
        'loudspeaker_model',            'Parametric', ...
        'loudspeaker_spacing',           0.01, ...
        'speaker_array_type',            'line'};
end

parametric_speaker = Parametric_Synthesis.parametric_soundfield;
parametric_speaker.P1 = db2mag( 100 ); % 100dB amplitude parametric array loudspeaker
parametric_speaker.P2 = db2mag( 100 ); % 100dB secondary amplitude

setup = Speaker_Setup.createSetup({...
    layout{:}, ...
    array{:}, ...
    'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
    'reproduction_radius',          1.0, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 0, ...
    'unattended_weight',            0.05, ...
    'brightzone_source_dist',       x_, ...
    'brightzone_radius',            0.3, ...
    'brightzone_pos_distance',      0.6, ...
    'quietzone_radius',             0.3, ...
    'quietzone_pos_distance',       0.6, ...
    'maximum_frequency',            8000, ...
    'angleof_loudspeakerarc',       180, ...
    'angleof_loudspeakerarrcentre',	180, ...
    'loudspeaker_object',           parametric_speaker});



Number_of_Frequencies = 512;

[~, Frequencies] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( ...
    Number_of_Frequencies, 0, 0, ...
    150, ...
    8000,'lin');

if spkr == 1
    Number_of_Weights = 32;
    Weights           = [0     logspace(log10( 1e-2), log10(  1e4), Number_of_Weights - 1) ];
elseif spkr == 2
    Number_of_Weights = 1;
    Weights           = 1;
end


%% Path to save the database (Look-Up Table)
DB_fullpath = Soundfield_Database.getDatabasePath( ...
    setup, ...
    [num2str(Number_of_Frequencies) 'f_' num2str(Number_of_Weights) 'w'], ...
    'Z:\', ...
    'new4');

DBpath = fileparts(DB_fullpath);

%% Results

%Samples
if spkr == 1
    Bright_Sample__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
    Quiet_Sample__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));
elseif spkr == 2
    Bright_Sample__Weight_Vs_Frequency = cell(length(Weights), length(Frequencies));
    Quiet_Sample__Weight_Vs_Frequency  = cell(length(Weights), length(Frequencies));
end

%Contrast
Acoustic_Contrast__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));

%Loudspeaker Weights
Loudspeaker_Weights__Weight_Vs_Frequency  = cell(length(Weights), length(Frequencies));


%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

fprintf('\n====== Building Look-Up Table (Frequency Weights) ======\n');
fprintf(Speaker_Setup.printSetupDetails(setup));
fprintf(['No. of Frequencies:\t' num2str(Number_of_Frequencies) '\n']);
fprintf(['No. of Weights:\t\t' num2str(Number_of_Weights) '\n\n']);
fprintf('\tCompletion: ');n=0;h=[];
a=1;
for w = 1:length(Weights)
    for f = 1:length(Frequencies)
        
        parsetup = setup;
        parsetup.Multizone_Soundfield.Quiet_Zone = ...
            parsetup.Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, Frequencies( f ), 'suppress_output');
        parsetup.Multizone_Soundfield.Bright_Zone = ...
            parsetup.Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, Frequencies( f ), 'suppress_output');
        
        parsetup.Multizone_Soundfield.QuietZ_Weight = Weights( w );
        parsetup.Multizone_Soundfield = parsetup.Multizone_Soundfield.setN( -1 ); %Auto set
        parsetup.Multizone_Soundfield = parsetup.Multizone_Soundfield.createSoundfield(DebugMode);
        
        parsetup = parsetup.calc_Loudspeaker_Weights();
        parsetup = parsetup.reproduceSoundfield('SAMPLES_ONLY');
        
        if spkr == 1
            Bright_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Bright_Sample;
            Quiet_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Quiet_Sample;
        elseif spkr == 2
            Bright_Sample__Weight_Vs_Frequency{ w, f } = parsetup.Bright_Samples;
            Quiet_Sample__Weight_Vs_Frequency{ w, f } = parsetup.Quiet_Samples;
        end
        Acoustic_Contrast__Weight_Vs_Frequency( w, f ) = parsetup.Acoustic_Contrast;
        Loudspeaker_Weights__Weight_Vs_Frequency{ w, f } = parsetup.Loudspeaker_Weights;
        
    end
    
    [n,h] = Tools.showTimeToCompletion(w/length(Weights), n, h);
    
end

Loudspeaker_Locations = setup.Loudspeaker_Locations;


%%

if ~exist(DBpath,'dir'), mkdir(DBpath);end;
save(DB_fullpath, ...
    'Frequencies', ...
    'Weights', ...
    'Bright_Sample__Weight_Vs_Frequency', ...
    'Quiet_Sample__Weight_Vs_Frequency', ...
    'Acoustic_Contrast__Weight_Vs_Frequency', ...
    'Loudspeaker_Weights__Weight_Vs_Frequency', ...
    'Loudspeaker_Locations');


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

end