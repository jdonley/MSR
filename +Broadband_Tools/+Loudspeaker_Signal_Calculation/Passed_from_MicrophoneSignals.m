function Passed_from_MicrophoneSignals( SYS )
% Generates spatially weighted loudspeaker signals passed from a microphone array.
% 
% Syntax:	Passed_from_MicrophoneSignals( SYS )
% 
% Inputs: 
% 	SYS - The Soundfield Reproduction system object
% 
% Example: 
% 	Passed_from_MicrophoneSignals(...
%           Current_Systems.loadCurrentSRsystem)
% 
% See also: getPredictedLoudspeakerSignals, loadCurrentSRsystem, createSetup

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 31 August 2017
% Version: 0.1 (31 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Determine Loudspeaker Signals
[Loudspeaker_Signals, Original_] = ...
    Broadband_Tools.Loudspeaker_Signal_Calculation.getPredictedLoudspeakerSignals( ...
    SYS );

%% Setup Variables

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, system_info.LUT_resolution, system_info.Drive, 'new');


%% Read input signal

if contains(lower(signal_info.method), 'cancel')
    %%% TODO
    % Fix this magnitude adjustment which is dependent on the transfer
    % function for point sources. Originates somewhere in
    % loudspeaker_setup.m or earlier. (Could be due to 2D ATFs and 3D RIRs)
    if setup.Dimensionality == 2
        magnitudeADJ = 0.8;
    else
        magnitudeADJ = 1.0;
    end
    Input_Signal = conv(-1,Input_Signal) * magnitudeADJ;
end


%% Normalise Loudspeaker Signals
if isfield(signal_info,'clipLevelAdjust') && isfield(signal_info,'inputSignalNorm') && signal_info.inputSignalNorm
    % Prevent clipping upon save
    Loudspeaker_Signals = Loudspeaker_Signals .* db2mag(signal_info.clipLevelAdjust);
end
Original_ = Original_ ./ max(abs(Original_(:)));


%% Once we have the speaker signals we should save them for later use as .wav files

Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    setup.Loudspeaker_Count, ...
    Original_, ...
    signal_info.input_filename, ...
    Input_file_ext, ...
    signal_info.Fs );



end