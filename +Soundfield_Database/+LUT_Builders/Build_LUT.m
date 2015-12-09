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
speech_layout = {'brightzone_pos_angle',        90, ...
                 'quietzone_pos_angle',         -90, ...
                 'brightzone_source_angle',     0};
masker_layout = {'brightzone_pos_angle',        -90, ...
                 'quietzone_pos_angle',         90, ...
                 'brightzone_source_angle',     -45};
             
setup = Speaker_Setup.createSetup({...
            masker_layout{:}, ...
            'resolution',                  100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
            'reproduction_radius',         1.0, ...
            'bright_weight',               1.0, ...
            'unattended_weight',           0.05, ...
            'brightzone_radius',           0.3, ...
            'brightzone_pos_distance',     0.6, ...
            'quietzone_radius',            0.3, ...
            'quietzone_pos_distance',      0.6, ...
            'numberof_loudspeakers',       24, ...
            'loudspeaker_radius',          1.5, ...
            'maximum_frequency',           8000, ...
            'angleto_firstloudspeaker',    90, ...
            'angleof_loudspeakerarc',      180, ...
            'loudspeaker_model',            'Genelec 8010A', ...
            'angleof_loudspeakerarrcentre',	180, ...
            'loudspeaker_spacing',          [], ...
            'speaker_array_type',           'circle'    });
    
                            
                            
Number_of_Frequencies = 512;
Number_of_Weights = 32;

[~, Frequencies] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( ...
                                     Number_of_Frequencies, 0, 0, ...
                                     150, ...
                                     8000,'lin');

Weights              = [0     logspace(log10( 1e-2), log10(  1e4), Number_of_Weights - 1) ];
%Weights              = [0     1e4];


%% Path to save the database (Look-Up Table)
DB_fullpath = Soundfield_Database.getDatabasePath( ...
    setup, ...
    [num2str(Number_of_Frequencies) 'f_' num2str(Number_of_Weights) 'w'], ...
    'Z:\', ...
    'new3');

DBpath = fileparts(DB_fullpath);

%% Results

%Samples
Bright_Sample__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Quiet_Sample__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));

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
    parfor f = 1:length(Frequencies)
        
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
                
        Bright_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Bright_Sample;
        Quiet_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Quiet_Sample;
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