%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison of fields from:
%  JIN, W., KLEIJN, W. B. & VIRETTE, D. Multizone soundfield reproduction 
%  using orthogonal basis expansion.  Acoustics, Speech and Signal Processing 
%  (ICASSP), 2013 IEEE International Conference on, 2013. IEEE, 311-315.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;
tic; %Start timing this script

%% Setup
resolution = 150; % Resolution of calculated fields (samples per metre)
f=2000;
Radius_Bright = 0.3;
Radius_Quiet = 0.3;
d = 0.6;
phi_d = 60;
phi_1 = 225;
phi_2 = 45;


%% Bright Zone
    bright = spatial_zone(f, Radius_Bright, 'pw', 1.0, phi_d);

%% Quiet Zone
    quiet = spatial_zone(f, Radius_Quiet, 'quiet');

%% Desired Soundfield
    Desired = global_soundfield;
    Desired = Desired.addSpatialZone(bright, d, phi_1);
    Desired = Desired.addSpatialZone(quiet,  d, phi_2);
    
%% Planewave Field
    Planewave = global_soundfield;
    Planewave = Planewave.addSpatialZone(bright, d, phi_1);

%% Build Fields
    Desired.res = resolution; Planewave.res = resolution; % Set resolution
    Desired   = Desired.createGlobalSoundfield;
    Planewave = Planewave.createGlobalSoundfield;
    
%% Plot Desired
    Desired.plotSoundfield;
    
%% Plot Planewave
    Planewave.plotSoundfield;
    
%% End    
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script    
 
    