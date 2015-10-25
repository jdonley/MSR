%% Step 0
% Generate soundfield databases
%Soundfield_Database.LUT_Builders.<some_database>

%% Step 1
% Generate Loudspeaker Signals from soundfield database
Broadband_Tools.BatchJob_Script__Generate_Loudspeaker_Signals;

%% Step 2
% Apply RIR's and save "recording"
Room_Acoustics.BatchJob_Script__Reverberant_MSR;

%% Step 3
% Load "recording", evaluate using different measures and save
Room_Acoustics.BatchJob_Script__Reverberant_MSR_Analysis;

%% Step 4
% Load results from evaluation and plot
Results.Plot_Reverb_PESQ_Results; % Includes Intelligiblity if set correctly