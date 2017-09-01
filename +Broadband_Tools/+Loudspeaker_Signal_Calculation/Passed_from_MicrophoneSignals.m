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
[Loudspeaker_Signals, Original] = ...
    Broadband_Tools.Loudspeaker_Signal_Calculation.getPredictedLoudspeakerSignals( ...
    SYS );

%% Setup Variables
lSYS = SYS; % Loudspeaker System (L-SYS)
sysI = strcmpi({lSYS.Room_Setup.SystemType},'transmit');
lSYS.Room_Setup = lSYS.Room_Setup( sysI );
lSYS.Main_Setup = lSYS.Main_Setup( sysI );
[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( ...
    lSYS.Main_Setup, ...
    SYS.signal_info, ...
    SYS.system_info.LUT_resolution, ...
    SYS.system_info.Drive );

mSYS = SYS; %  Microphone System (M-SYS)
sysI = strcmpi({mSYS.Room_Setup.SystemType},'receive');
mSYS.Room_Setup = mSYS.Room_Setup( sysI );
mSYS.Main_Setup = mSYS.Main_Setup( sysI );
mpath = Broadband_Tools.getMicrophoneSignalPath( mSYS );
mfiles = Tools.keepFilesFromFolder( Tools.getAllFiles(mpath), mSYS.signal_info.speech_filepath);
[~,~,Input_file_ext] = fileparts(mfiles{contains(lower(mfiles),'original')});

%% Read input signal

if contains(lower(SYS.signal_info.method), 'cancel')
    %%% TODO
    % Fix this magnitude adjustment which is dependent on the transfer
    % function for point sources. Originates somewhere in
    % loudspeaker_setup.m or earlier. (Could be due to 2D ATFs and 3D RIRs)
    if lSYS.Main_Setup.Dimensionality == 2
        magnitudeADJ = 0.8;
    else
        magnitudeADJ = 1.0;
    end
    Loudspeaker_Signals = -Loudspeaker_Signals * magnitudeADJ;
end


%% Normalise Loudspeaker Signals
if isfield(SYS.signal_info,'clipLevelAdjust') ...
        && isfield(SYS.signal_info,'inputSignalNorm') ...
        && SYS.signal_info.inputSignalNorm
    % Prevent clipping upon save
    Loudspeaker_Signals = Loudspeaker_Signals .* db2mag(SYS.signal_info.clipLevelAdjust);
end
Original = Original ./ max(abs(Original(:)));


%% Once we have the speaker signals we should save them for later use as .wav files

Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    lSYS.Main_Setup.Loudspeaker_Count, ...
    Original, ...
    SYS.signal_info.input_filename, ...
    Input_file_ext, ...
    SYS.signal_info.Fs );



end