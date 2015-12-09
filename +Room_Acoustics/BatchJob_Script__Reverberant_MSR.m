delete(gcp);clear;close all;fclose all;clc;clear classes;

%%
%rooms = [4];
rooms = [ 1 ];

%setups = [1 2 3 4 5 6 0];
setups = [ 7.7 ];

for room = rooms
    
    Noise_Mask_Levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20];% 25 30 35 40];
    loudspeaker_layout = { ...
        'numberof_loudspeakers',        24, ...
        'loudspeaker_radius',           1.5, ...
        'loudspeaker_model',            'Genelec 8010A', ...
        'angleof_loudspeakerarrcentre', 180, ...
        'loudspeaker_spacing',          [], ...
        'speaker_array_type',           'circle' };
    speech_layout = {};
    masker_layout = {};
    
    Room_Setup = Room_Acoustics.Room;
    Room_Setup.NoReceivers = 32;
    
    %%
    if room == 1
        % % ROOM 1
        % % Anechoic
        Room_Setup.Room_Size = [10 10 10]; %Anechoic
        Room_Setup.Wall_Absorb_Coeff = 1.0;
        
    elseif room == 2
        % % ROOM 2
        % % Small Open Plan Office
        Room_Setup.Room_Size = [4 9 3];   %Small Open Plan Office
        Room_Setup.Wall_Absorb_Coeff = 0.3;
        
    elseif room == 3
        % % ROOM 3
        % % Medium Open Plan Office
        Room_Setup.Room_Size = [8 10 3];   %Medium Open Plan Office
        Room_Setup.Wall_Absorb_Coeff = 0.3;
        
    elseif room == 4
        % % ROOM 4
        % % Cafe / Restaurant
        Room_Setup.Room_Size = [9 14 3];   %Cafe/Restaurant
        Room_Setup.Wall_Absorb_Coeff = 0.3;
    end
    
    %%
    for scheme = setups
        if scheme == 0
            % % Setup and Privacy Scheme 0 (Benchmark)
            %Zone_Weights = 0;
            mask_type = 'NoMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     0};
            Noise_Mask_Levels = [];
        elseif scheme == 1
            % % Setup and Privacy Scheme 1
            mask_type = 'FlatMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     0};
        elseif scheme == 2
            % % Setup and Privacy Scheme 2
            Masker_Type = 'ZoneWeightMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     0};
        elseif scheme == 3
            % % Setup and Privacy Scheme 3
            mask_type = 'FlatMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     15};
        elseif scheme == 4
            % % Setup and Privacy Scheme 4
            mask_type = 'ZoneWeightMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     15};
        elseif scheme == 5
            % % Setup and Privacy Scheme 5
            mask_type = 'FlatMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     90};
        elseif scheme == 6
            % % Setup and Privacy Scheme 6
            mask_type = 'ZoneWeightMask';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     90};
            
            %Test Schemes
        elseif scheme == 7
            % % Setup and Privacy Scheme 7
            mask_type = 'ZoneWeightMaskAliasCtrl';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     15};
        elseif scheme == 7.5
            % % Setup and Privacy Scheme 7.5
            mask_type = 'ZoneWeightMaskAliasCtrl';
            speech_layout = {'brightzone_pos_angle',        90, ...
                'quietzone_pos_angle',         -90, ...
                'brightzone_source_angle',     0};
        elseif scheme == 7.6
            % % Setup and Privacy Scheme 7.6
            mask_type = 'ZoneWeightMaskAliasCtrl';
            speech_layout = {'brightzone_pos_angle',        90, ...
                'quietzone_pos_angle',         -90, ...
                'brightzone_source_angle',     0};
        elseif scheme == 7.7
            % % Setup and Privacy Scheme 7.7
            mask_type = 'ZoneWeightMaskAliasCtrl';
            speech_layout = {'brightzone_pos_angle',        90, ...
                'quietzone_pos_angle',         -90, ...
                'brightzone_source_angle',     0};
        elseif scheme == 8
            % % Setup and Privacy Scheme 8
            mask_type = 'ZoneWeightMaskAliasCtrlOffsetNoise';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     15};
        elseif scheme == 9
            % % Setup and Privacy Scheme 9
            mask_type = 'ZoneWeightMaskAliasCtrlStereoNoise';
            speech_layout = {'brightzone_pos_angle',        180, ...
                'quietzone_pos_angle',         0, ...
                'brightzone_source_angle',     15};
        end
        
        Speech_Setup = Speaker_Setup.createSetup({ speech_layout{:}, loudspeaker_layout{:}});
        
        Room_Acoustics.Apply_RIRs.Reverberant_MSR_batchfunc( Room_Setup, mask_type, Speech_Setup, Noise_Mask_Levels);
        
    end
end