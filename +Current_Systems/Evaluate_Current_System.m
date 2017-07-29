
%                       COMPLETE EVALUATION SCRIPT
%--------------------------------------------------------------------------
% Make sure the current system can be loaded correctly prior to simulation 
%==========================================================================
clc;clear;clear classes;close all;fclose all;delete(gcp('nocreate'));

%% Build Look-Up Tables
Soundfield_Database.LUT_Builders.Build_LUT

%% Build Room Impulse Responses
Room_Acoustics.Generate_RIR_Database


%% Generate Loudspeaker Signals
Broadband_Tools.Generate_Loudspeaker_Signals


%% Obtain Calibration Filters
Speaker_Setup.Calibration.Perform_Full_Calibration

%% Calibrate Loudspeaker Signals
Speaker_Setup.Calibration.Calibrate_MSR_Loudspeaker_Signals

%% Realworld Recording
Hardware_Control.Play_and_Rec_System


%% Simulate Recording 
Room_Acoustics.Reverberant_MSR

%% Measure Recorded Signal Characteristics
Room_Acoustics.Reverberant_MSR_Analysis

%% Plot Results
Results.Plot_Results
