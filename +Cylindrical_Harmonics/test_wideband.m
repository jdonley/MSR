clc;
clear;
close all;

%% Setup
    zone = [];
    monofield = [];
    soundfield = [];


%% Zone 1, Frequency 1
    % Sound in Zone: Frequency, Radius, Angle, Source Type, (Source Distance)
    zone = [zone spatial_zone(4000, 0.3, 0, 'pw')]; 

%% Zone 1, Frequency 2
    % Sound in Zone: Frequency, Radius, Angle, Source Type, (Source Distance)
    zone = [zone spatial_zone(700, 0.3, 0, 'pw')];

%% Zone 2, Frequency 1
    % Sound in Zone: Frequency, Radius, Angle, Source Type, (Source Distance)
    zone = [zone spatial_zone(2000, 0.3, 10, 'pw')]; 

%% Zone 2, Frequency 2
    % Sound in Zone: Frequency, Radius, Angle, Source Type, (Source Distance)
    zone = [zone spatial_zone(300, 0.3, 10, 'pw')];


%% Frequency 1
    monofield = [monofield global_soundfield];
    % Location of Zone: Distance and Angle from Origin
    monofield(1) = monofield(1).addSpatialZone(zone(1), 1.5, 0);

%% Frequency 2
    monofield = [monofield global_soundfield];
    % Location of Zone: Distance and Angle from Origin
    monofield(2) = monofield(2).addSpatialZone(zone(2), 1.5, 0);

%% Frequency 3
    monofield = [monofield global_soundfield];
    % Location of Zone: Distance and Angle from Origin
    monofield(3) = monofield(3).addSpatialZone(zone(3), 1.5, 120);
    
%% Frequency 4
    monofield = [monofield global_soundfield];
    % Location of Zone: Distance and Angle from Origin
    monofield(4) = monofield(4).addSpatialZone(zone(4), 1.5, 120);


%% Wideband field
    soundfield = wideband_soundfield;
    %Add each monofield to the wideband field
    for f = 1:length(monofield)
        soundfield = soundfield.addGlobalSoundfield(monofield(f));
    end
    
    %Create wideband field
    soundfield = soundfield.setResolution(50); % Resolution is in samples per metre
    soundfield = soundfield.createWidebandSoundfield();


%% Display Sound Field
    soundfield.plotSoundfield;