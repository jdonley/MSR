delete(gcp);clear;close all;fclose all;clc;

%%
%rooms = [4];
rooms = [ 1 ];

%setups = [1 2 3 4 5 6 0];
setups = [ 9 ];

for room = rooms
    
Noise_Mask_Levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20];% 25 30 35 40];
    
    
%%
    if room == 1
        % % ROOM 1
        % % Anechoic
        Room_Size = [10 10 10]; %Anechoic
        Wall_Absorption_Coeff = 1.0;
        
    elseif room == 2
        % % ROOM 2
        % % Small Open Plan Office
        Room_Size = [4 9 3];   %Small Open Plan Office
        Wall_Absorption_Coeff = 0.3;
        
    elseif room == 3
        % % ROOM 3
        % % Medium Open Plan Office
        Room_Size = [8 10 3];   %Medium Open Plan Office
        Wall_Absorption_Coeff = 0.3;
        
    elseif room == 4
        % % ROOM 4
        % % Cafe / Restaurant
        Room_Size = [9 14 3];   %Cafe/Restaurant
        Wall_Absorption_Coeff = 0.3;
    end    
    
    %%    
    for scheme = setups
        if scheme == 0
            % % Setup and Privacy Scheme 0 (Benchmark)
            %Zone_Weights = 0;
            mask_type = 'NoMask';
            pw_angle = 0;
            Noise_Mask_Levels = [];
        elseif scheme == 1
            % % Setup and Privacy Scheme 1
            mask_type = 'FlatMask';
            pw_angle = 0;
        elseif scheme == 2
            % % Setup and Privacy Scheme 2
            mask_type = 'ZoneWeightMask';
            pw_angle = 0;
        elseif scheme == 3
            % % Setup and Privacy Scheme 3
            mask_type = 'FlatMask';
            pw_angle = 15;
        elseif scheme == 4
            % % Setup and Privacy Scheme 4
            mask_type = 'ZoneWeightMask';
            pw_angle = 15;
        elseif scheme == 5
            % % Setup and Privacy Scheme 5
            mask_type = 'FlatMask';
            pw_angle = 90;
        elseif scheme == 6
            % % Setup and Privacy Scheme 6
            mask_type = 'ZoneWeightMask';
            pw_angle = 90;
        
            %Test Schemes
        elseif scheme == 7
            % % Setup and Privacy Scheme 7
            mask_type = 'ZoneWeightMaskAliasCtrl';
            pw_angle = 15;
        elseif scheme == 8
            % % Setup and Privacy Scheme 7
            mask_type = 'ZoneWeightMaskAliasCtrlOffsetNoise';
            pw_angle = 15;
        elseif scheme == 9
            % % Setup and Privacy Scheme 7
            mask_type = 'ZoneWeightMaskAliasCtrlStereoNoise';
            pw_angle = 15;
        end
        
        Room_Acoustics.Apply_RIRs.Reverberant_MSR_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, pw_angle, Noise_Mask_Levels);
        
    end
end