clear;close all;fclose all;clc;

Zone_Weights = 1e4;
Noise_Mask_Levels = [-40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30 35 40];


%%
setups = [1 2 3 4 5 6];

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
    end
    
    Broadband_Tools.Generate_Loudspeaker_Signals_batchfunc(Zone_Weights, Noise_Mask_Levels, Masker_Type, Planewave_Angle);
    
end