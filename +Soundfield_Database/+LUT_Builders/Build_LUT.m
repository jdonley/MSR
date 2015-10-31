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

speech_layout = {'brightzone_pos_angle',        180, ...
                 'quietzone_pos_angle',         0, ...
                 'brightzone_source_angle',     15};
masker_layout = {'brightzone_pos_angle',        0, ...
                 'quietzone_pos_angle',         180, ...
                 'brightzone_source_angle',     0};
             
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
            'loudspeaker_spacing',          0.01    });
    
                            
                            

Number_of_Frequencies = 512;
Number_of_Weights = 32;

[NumberOfPlanewaves, Frequencies] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( ...
                                     Number_of_Frequencies, ...
                                     28, ...
                                     300, ...
                                     150, ...
                                     8000,'lin');

Weights              = [0     logspace(log10( 1e-2), log10(  1e4), Number_of_Weights - 1) ];
%Weights              = [0     1e4];

SPL = 94; % 94dB is 1Pa (RMS) when reference to 20uPa (RMS)


% Path to save the database (Look-Up Table)
Database_Path = ['Z:\+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];


%% Results
Contrast__Weight_Vs_Frequency     = zeros(length(Weights), length(Frequencies));
Error_Bright__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Error_Quiet__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));
Quiet_SPL__Weight_Vs_Frequency    = zeros(length(Weights), length(Frequencies));

%Samples
Bright_Sample__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Quiet_Sample__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));

%Loudspeaker Weights
Loudspeaker_Weights__Weight_Vs_Frequency  = cell(length(Weights), length(Frequencies));


%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

fprintf('\n====== Building Look-Up Table (Frequency Weights) ======\n\n');
fprintf('\tCompletion: ');n=0;
a=1;
for w = 1:length(Weights)
    parfor f = 1:length(Frequencies)
        
        parsetup = setup;
        parsetup.Multizone_Soundfield.Quiet_Zone = ...
            parsetup.Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, Frequencies( f ), 'suppress_output');
        parsetup.Multizone_Soundfield.Bright_Zone = ...
            parsetup.Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, Frequencies( f ), 'suppress_output');
        
        parsetup.Multizone_Soundfield.QuietZ_Weight = Weights( w );
        parsetup.Multizone_Soundfield = parsetup.Multizone_Soundfield.setN( NumberOfPlanewaves( f ) );
        parsetup.Multizone_Soundfield = parsetup.Multizone_Soundfield.createSoundfield(DebugMode);
                       
        parsetup = parsetup.calc_Loudspeaker_Weights();
        parsetup = parsetup.reproduceSoundfield('SAMPLES_ONLY');
                
        Bright_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Bright_Sample;
        Quiet_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Quiet_Sample;                
        Loudspeaker_Weights__Weight_Vs_Frequency{ w, f } = parsetup.Loudspeaker_Weights;                

    end
    
    n = Tools.showTimeToCompletion(w/length(Weights), n);
    
end

Loudspeaker_Locations = setup.Loudspeaker_Locations;


%%
if ~exist(Database_Path,'dir'), mkdir(Database_Path);end;
save([Database_Path 'LUT_Weight_vs_Frequency_' ...
      num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' ...
      num2str(Number_of_Frequencies) 'f_' ...
      num2str(Number_of_Weights) 'w.mat'], ...
     'Frequencies', ...
     'Weights', ...
     'Error_Quiet__Weight_Vs_Frequency', ...
     'Error_Bright__Weight_Vs_Frequency', ...
     'Contrast__Weight_Vs_Frequency', ...
     'Quiet_SPL__Weight_Vs_Frequency', ...
     'Bright_Sample__Weight_Vs_Frequency', ...
     'Quiet_Sample__Weight_Vs_Frequency', ...
     'Loudspeaker_Weights__Weight_Vs_Frequency', ...
     'Loudspeaker_Locations');


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script