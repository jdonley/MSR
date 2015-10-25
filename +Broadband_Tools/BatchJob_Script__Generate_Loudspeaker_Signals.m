clear;close all;fclose all;clc;

Zone_Weights = 1e4;
Noise_Mask_Levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20];% 25 30 35 40];
Leakage_Angle = 0;

%%
setups = [7];
%setups = [0];

for scheme = setups
    if scheme == 0
        % % Setup and Privacy Scheme 0 (Benchmark)
        Zone_Weights = 0;
        Masker_Type = 'NoMask';        
        Planewave_Angle = 0;
        Noise_Mask_Levels = -inf;
    elseif scheme == 1
        % % Setup and Privacy Scheme 1
        Masker_Type = 'FlatMask';
        Planewave_Angle = 0;
    elseif scheme == 2
        % % Setup and Privacy Scheme 2
        Masker_Type = 'ZoneWeightMask';
        Planewave_Angle = 0;
    elseif scheme == 3
        % % Setup and Privacy Scheme 3
        Masker_Type = 'FlatMask';
        Planewave_Angle = 15;
    elseif scheme == 4
        % % Setup and Privacy Scheme 4
        Masker_Type = 'ZoneWeightMask';
        Planewave_Angle = 15;
    elseif scheme == 5
        % % Setup and Privacy Scheme 5
        Masker_Type = 'FlatMask';
        Planewave_Angle = 90;
    elseif scheme == 6
        % % Setup and Privacy Scheme 6
        Masker_Type = 'ZoneWeightMask';
        Planewave_Angle = 90;
        
    %Test Schemes    
    elseif scheme == 7
        % % Setup and Privacy Scheme 7
        Masker_Type = 'ZoneWeightMask_AliasCtrl';
        Planewave_Angle = 15;
    elseif scheme == 8
        % % Setup and Privacy Scheme 8
        Masker_Type = 'ZoneWeightMask_AliasCtrl_OffsetNoise';
        Planewave_Angle = 15;
        Leakage_Angle = 15;
    elseif scheme == 9
        % % Setup and Privacy Scheme 9
        Masker_Type = 'ZoneWeightMask_AliasCtrl_StereoNoise';
        Planewave_Angle = 15;
        Leakage_Angle = 30;
    end
    
    Broadband_Tools.Generate_Loudspeaker_Signals_batchfunc(Zone_Weights, Noise_Mask_Levels, Masker_Type, Planewave_Angle, Leakage_Angle);
    
end