clear;close all;fclose all;clc;

Zone_Weights = 1e4;
Noise_Mask_Levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20];% 25 30 35 40];
Leakage_Angle = 0;
loudspeaker_layout = {'numberof_loudspeakers',        24, ...
                      'loudspeaker_radius',           1.5, ...
                      'loudspeaker_model',            'Genelec 8010A', ...
                      'angleof_loudspeakerarrcentre', 180, ...
                      'loudspeaker_spacing',          [], ...
                      'speaker_array_type',           'circle' };
speech_layout = {};
masker_layout = {};

%%
setups = [7.7];
%setups = [0];

for scheme = setups
    if scheme == 0
        % % Setup and Privacy Scheme 0 (Benchmark)
        Zone_Weights = 0;
        Masker_Type = 'NoMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     0};
        Noise_Mask_Levels = -inf;
    elseif scheme == 1
        % % Setup and Privacy Scheme 1
        Masker_Type = 'FlatMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     0};
    elseif scheme == 2
        % % Setup and Privacy Scheme 2
        Masker_Type = 'ZoneWeightMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     0};
        masker_layout = {'brightzone_pos_angle',        0, ...
                         'quietzone_pos_angle',         180, ...
                         'brightzone_source_angle',     0};
    elseif scheme == 3
        % % Setup and Privacy Scheme 3
        Masker_Type = 'FlatMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     15};
    elseif scheme == 4
        % % Setup and Privacy Scheme 4
        Masker_Type = 'ZoneWeightMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     15};
        masker_layout = {'brightzone_pos_angle',        0, ...
                         'quietzone_pos_angle',         180, ...
                         'brightzone_source_angle',     0};
    elseif scheme == 5
        % % Setup and Privacy Scheme 5
        Masker_Type = 'FlatMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     90};
    elseif scheme == 6
        % % Setup and Privacy Scheme 6
        Masker_Type = 'ZoneWeightMask';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     90};
        masker_layout = {'brightzone_pos_angle',        0, ...
                         'quietzone_pos_angle',         180, ...
                         'brightzone_source_angle',     0};
        
    %Test Schemes    
    elseif scheme == 7
        % % Setup and Privacy Scheme 7
        Masker_Type = 'ZoneWeightMask_AliasCtrl';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     15};
        masker_layout = {'brightzone_pos_angle',        0, ...
                         'quietzone_pos_angle',         180, ...
                         'brightzone_source_angle',     0};
    elseif scheme == 7.5
        % % Setup and Privacy Scheme 7.5
        Masker_Type = 'ZoneWeightMask_AliasCtrl';
        speech_layout = {'brightzone_pos_angle',        90, ...
                         'quietzone_pos_angle',         -90, ...
                         'brightzone_source_angle',     0};
        masker_layout = {'brightzone_pos_angle',        -90, ...
                         'quietzone_pos_angle',         90, ...
                         'brightzone_source_angle',     -60};
    elseif scheme == 7.6
        % % Setup and Privacy Scheme 7.6
        Masker_Type = 'ZoneWeightMask_AliasCtrl';
        speech_layout = {'brightzone_pos_angle',        90, ...
                         'quietzone_pos_angle',         -90, ...
                         'brightzone_source_angle',     0};
        masker_layout = {'brightzone_pos_angle',        -90, ...
                         'quietzone_pos_angle',         90, ...
                         'brightzone_source_angle',     0};
    elseif scheme == 7.7
        % % Setup and Privacy Scheme 7.7
        Masker_Type = 'ZoneWeightMask_AliasCtrl';
        speech_layout = {'brightzone_pos_angle',        90, ...
                         'quietzone_pos_angle',         -90, ...
                         'brightzone_source_angle',     0};
        masker_layout = {'brightzone_pos_angle',        -90, ...
                         'quietzone_pos_angle',         90, ...
                         'brightzone_source_angle',     -45};
    elseif scheme == 8
        % % Setup and Privacy Scheme 8
        Masker_Type = 'ZoneWeightMask_AliasCtrl_OffsetNoise';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     15};
        masker_layout = {'brightzone_pos_angle',        0, ...
                         'quietzone_pos_angle',         180, ...
                         'brightzone_source_angle',     15};
    elseif scheme == 9
        % % Setup and Privacy Scheme 9
        Masker_Type = 'ZoneWeightMask_AliasCtrl_StereoNoise';
        speech_layout = {'brightzone_pos_angle',        180, ...
                         'quietzone_pos_angle',         0, ...
                         'brightzone_source_angle',     15};
        masker_layout = {'brightzone_pos_angle',        0, ...
                         'quietzone_pos_angle',         180, ...
                         'brightzone_source_angle',     30};
    end
    
    Speech_Setup = Speaker_Setup.createSetup({ speech_layout{:}, loudspeaker_layout{:}});    
    Masker_Setup = Speaker_Setup.createSetup({ masker_layout{:}, loudspeaker_layout{:}});
    
    Broadband_Tools.Generate_Loudspeaker_Signals_batchfunc(Zone_Weights, Noise_Mask_Levels, Masker_Type, Speech_Setup, Masker_Setup);
    
end