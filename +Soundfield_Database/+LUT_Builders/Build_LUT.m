function Build_LUT( SYS )

%clear;
close all;
tic; %Start timing this script
delete(gcp('nocreate'));
current_pool = parpool; %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))



%% Settings
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end



DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
% You can suppress the output when executing a complete soundfield construction with 'suppress_output'



%%
Setups = [];
if ~isempty(SYS.signal_info.methods_list_clean)
    Setups = [Setups(:); SYS.Main_Setup(:)];
end
if ~isempty(SYS.signal_info.methods_list_masker)
    Setups = [Setups(:); SYS.Masker_Setup(:)];
end

DBsetups = 1:length(Setups);
if isfield(SYS.system_info,'DB_indices')
    DBsetups(~reshape([SYS.system_info.DB_indices{:}],1,[]))=[];
end
for s = DBsetups
    Setup = Setups(s);
    
    [~, Frequencies] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( ...
        SYS.system_info.LUT_frequencies, 0, 0, ...
        SYS.signal_info.f_low, ...
        SYS.signal_info.f_high,'lin');
    
    single_LUT_weight = (SYS.system_info.LUT_weights == 1);
    if Setup.Loudspeaker_Count > 1
        Weights = [0,  logspace(log10( min(SYS.system_info.LUT_weight_range) ), ...
            log10( max(SYS.system_info.LUT_weight_range) ), ...
            SYS.system_info.LUT_weights - 1) ];
    elseif Setup.Loudspeaker_Count == 1 && single_LUT_weight
        Weights = 1;
    end
    
    
    %% Path to save the database (Look-Up Table)
    DB_fullpath = Soundfield_Database.getDatabasePath( ...
        Setup, ...
        [num2str(SYS.system_info.LUT_frequencies) 'f_' num2str(length(Weights)) 'w'], ...
        SYS.system_info.Drive);
    
    DBpath = fileparts(DB_fullpath);
    
    %% Results
    
    %Samples
    if Setup.Loudspeaker_Count > 1
        Bright_Sample__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
        Quiet_Sample__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));
    elseif Setup.Loudspeaker_Count == 1 && single_LUT_weight
        Bright_Sample__Weight_Vs_Frequency = cell(length(Weights), length(Frequencies));
        Quiet_Sample__Weight_Vs_Frequency  = cell(length(Weights), length(Frequencies));
    end
    
    %Contrast
    Acoustic_Contrast__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
    
    %Loudspeaker Weights
    Loudspeaker_Weights__Weight_Vs_Frequency  = cell(length(Weights), length(Frequencies));
    
    
    %% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)
    
    fprintf('\n====== Building Look-Up Table (Frequency Weights) ======\n');
    fprintf(Speaker_Setup.printSetupDetails(Setup));
    fprintf(['No. of Frequencies:\t' num2str(SYS.system_info.LUT_frequencies) '\n']);
    fprintf(['No. of Weights:\t\t' num2str(length(Weights)) '\n\n']);
    fprintf('\t Completion: ');n=0;
    a=1;
    for w = 1:length(Weights)
        currWeight = Weights( w ); % Current Weight to process
        parfor f = 1:length(Frequencies)
            
            parsetup = Setup;
            parsetup.Multizone_Soundfield.Quiet_Zone = ...
                parsetup.Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, Frequencies( f ), 'suppress_output');
            parsetup.Multizone_Soundfield.Bright_Zone = ...
                parsetup.Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, Frequencies( f ), 'suppress_output');
            
            parsetup.Multizone_Soundfield.QuietZ_Weight = currWeight;
            parsetup.Multizone_Soundfield = parsetup.Multizone_Soundfield.setN( -1 ); %Auto set
            parsetup.Multizone_Soundfield = parsetup.Multizone_Soundfield.createSoundfield(DebugMode);
            
            parsetup = parsetup.calc_Loudspeaker_Weights();
            parsetup = parsetup.reproduceSoundfield('SAMPLES_ONLY');
            
            if parsetup.Loudspeaker_Count > 1
                Bright_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Bright_Sample;
                Quiet_Sample__Weight_Vs_Frequency( w, f ) = parsetup.Quiet_Sample;
            elseif parsetup.Loudspeaker_Count == 1 && single_LUT_weight
                Bright_Sample__Weight_Vs_Frequency{ w, f } = parsetup.Bright_Samples;
                Quiet_Sample__Weight_Vs_Frequency{ w, f } = parsetup.Quiet_Samples;
            end
            Acoustic_Contrast__Weight_Vs_Frequency( w, f ) = parsetup.Acoustic_Contrast;
            Loudspeaker_Weights__Weight_Vs_Frequency{ w, f } = parsetup.Loudspeaker_Weights;
            
        end
        
        n = Tools.showTimeToCompletion(w/length(Weights), n);
        
    end
    
    Loudspeaker_Locations = Setup.Loudspeaker_Locations;
    
    
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