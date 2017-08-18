function Clean_from_RIR( Input_file, SYS )
% Generates microphone signal recordings for an array of microphones
% 
% Syntax:	Clean_from_RIR( Input_file, SYS_or_LUT_resolution, setup )
% 
% Inputs: 
% 	Input_file - The input audio file to generate loudspeaker signals for.
% 	SYS_or_LUT_resolution - Soundfield Reproduction system object
% 
% Example: 
% 	Clean_from_RIR('C:\someaudiofile.wav', ...
%       Current_Systems.loadCurrentSRsystem)
% 
% See also: loadCurrentSRsystem, createSetup

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 17 August 2017
% Version: 0.1 (17 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup Variables
if isempty(SYS.signal_info.method)
    SYS.signal_info.method = 'Clean';
end
room = SYS.Room_Setup;
signal_info = SYS.signal_info;
system_info = SYS.system_info;
setup = SYS.Main_Setup;
    
[~, Input_file_name, Input_file_ext] = fileparts( Input_file );
signal_info.input_filename = Input_file_name;



[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getMicrophoneSignalPath( SYS );



%% Read input signal
Input_Signal = audioread( Input_file );

if isfield(signal_info,'inputSignalNorm') && signal_info.inputSignalNorm
    % Normalise input signal
    Input_Signal = Input_Signal ./ rms(Input_Signal);
end

%% Generate Loudspeaker Signals
[Microphone_Signals, Original_] = Broadband_Tools.Microphone_Signal_Calculation.getMicrophoneSignals( ...
    Input_Signal, SYS );

%% Normalise Loudspeaker Signals
if isfield(signal_info,'clipLevelAdjust') && isfield(signal_info,'inputSignalNorm') && signal_info.inputSignalNorm
    % Prevent clipping upon save
    Microphone_Signals = Microphone_Signals .* db2mag(signal_info.clipLevelAdjust);
end
Original_ = Original_ ./ max(abs(Original_(:)));




%% Once we have the speaker signals we should save them for later use as .wav files

Broadband_Tools.Microphone_Signal_Calculation.saveMicrophoneSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Microphone_Signals, ...
    room.NoReceivers, ...
    Original_, ...
    signal_info.input_filename, ...
    Input_file_ext, ...
    signal_info.Fs );



end