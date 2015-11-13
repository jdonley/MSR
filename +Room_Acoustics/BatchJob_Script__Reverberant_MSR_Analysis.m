delete(gcp);clear;close all;fclose all;clc;

% % ROOM 1
% % Anechoic6
pesqNumber = 0;
Room_Size = [10 10 10]; %Anechoic
Wall_Absorption_Coeff = 1.0;

% % ROOM 2
% % Small Open Plan Office
% pesqNumber = 1;
% Room_Size = [4 9 3];   %Small Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 3
% % Medium Open Plan Office
% pesqNumber = 2;
% Room_Size = [8 10 3];   %Medium Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 4
% % Cafe / Restaurant
% pesqNumber = 4;
% Room_Size = [9 14 3];   %Cafe/Restaurant
% Wall_Absorption_Coeff = 0.3;

%%
%setups = [0];
setups = [ 7.6 ];
%% 
Noise_Mask_Levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20];% 25 30 35 40];
    loudspeaker_layout = {'numberof_loudspeakers',        24, ...
        'loudspeaker_radius',           1.5, ...
        'loudspeaker_model',            'Genelec 8010A', ...
        'angleof_loudspeakerarrcentre', 180, ...
        'loudspeaker_spacing',          0.01, ...
        'speaker_array_type',           'line'    };
    speech_layout = {};
    masker_layout = {};
    

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
        
    Room_Acoustics.Apply_RIRs.Reverberant_MSR_Analysis_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, Speech_Setup, Noise_Mask_Levels, pesqNumber);
    %Room_Acoustics.Apply_RIRs.Reverberant_MSR_Analysis_STI_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, pw_angle, Noise_Mask_Levels);
    
end