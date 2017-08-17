
%                       COMPLETE EVALUATION SCRIPT
%--------------------------------------------------------------------------
% Make sure the current system can be loaded correctly prior to simulation 
%==========================================================================
clc;clear;clear classes;close all;fclose all;delete(gcp('nocreate'));

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016-2017
% Date: 17 August 2017
% Version: 0.2 (17 August 2017)
% Version: 0.1 (14 June 2016)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Build Receiver and Reproduction Room Impulse Responses
Room_Acoustics.Generate_RIR_Database

%% Generate Microphone Signals
Broadband_Tools.Generate_Microphone_Signals


%% Build Reproduction Look-Up Tables
Soundfield_Database.LUT_Builders.Build_LUT

%% Generate Loudspeaker Signals
Broadband_Tools.Generate_Loudspeaker_Signals


%% Obtain Calibration Filters
Speaker_Setup.Calibration.Perform_Full_Calibration

%% Calibrate Loudspeaker Signals
Speaker_Setup.Calibration.Calibrate_MSR_Loudspeaker_Signals


%% Realworld Reproduction Recording
Hardware_Control.Play_and_Rec_System

%% Simulated Reproduction Recording 
Room_Acoustics.Reverberant_MSR


%% Measure Recorded Signal Characteristics
Room_Acoustics.Reverberant_MSR_Analysis

%% Plot Results
Results.Plot_Results
