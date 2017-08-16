function Clean_from_LUT( Input_file, SYS_or_LUT_resolution, weight, setup )
% Generates spatially weighted loudspeaker signals for the reproduction of zoned clean speech.
% 
% Syntax:	CLEAN_FROM_LUT( INPUT_FILE, SYS_OR_LUT_RESOLUTION, WEIGHT, SETUP )
% 
% Inputs: 
% 	Input_file - The input audio file to generate loudspeaker signals for.
% 	SYS_or_LUT_resolution - Either a Soundfield Reproduction system object
% 	or a string containing the LUT resolution in the format:
%   '<NFrequencies>f_<Nweights>w', where <NFrequencies> is the number of
%   frequency coefficients and <Nweights> is the number of different zone
%   weights.
% 	weight - The multizone soundfield weight
% 	setup - Loudspeaker setup object
% 
% Example: 
% 	Clean_from_LUT('C:\someaudiofile.wav', ...
%       Current_Systems.loadCurrentSRsystem)
% 
% See also: getMSRLoudspeakerSignals, loadCurrentSRsystem, createSetup

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 15 August 2015
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SYS_type = 'Current_Systems.SR_System';

%% Setup Variables
if nargin == 2
    if isa(SYS_or_LUT_resolution,SYS_type)
        SYS = SYS_or_LUT_resolution;        
        signal_info = SYS.signal_info;
        system_info = SYS.system_info;
        setup = SYS.Main_Setup;
    else
        error(['Second input argument must be of type: ' SYS_type]);
    end
elseif nargin == 4
    system_info.Drive = 'Z:\';
    system_info.LUT_resolution = SYS_or_LUT_resolution;
    
    signal_info.c = 343; % Speed of sound in metres/sec
    signal_info.Fs = 16000; % Sampling frequency
    signal_info.Nfft = 1024;% Number of fft components
    signal_info.overlap = 0.5;
    signal_info.f_low  = 150;  % Hz
    signal_info.f_high = 8000; % Hz
    signal_info.f_high_meas = 7000; % Hz %Maximum frequency with accurate response at given sampling rate
    signal_info.weight = weight;
    signal_info.method = 'NoMask';
end

% Clean speech means no noise levels
signal_info.L_noise_mask = -Inf; % dB
    
[~, Input_file_name, Input_file_ext] = fileparts( Input_file );
signal_info.input_filename = Input_file_name;



[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, system_info.LUT_resolution, system_info.Drive, 'new');



%% Read input signal
Input_Signal = audioread( Input_file );

if isfield(signal_info,'inputSignalNorm') && signal_info.inputSignalNorm
    % Normalise input signal
    Input_Signal = Input_Signal ./ rms(Input_Signal);
end

if ~isempty(strfind(lower(signal_info.method), 'cancel'))
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

%% Generate Loudspeaker Signals
[Loudspeaker_Signals, Original_] = Broadband_Tools.Loudspeaker_Signal_Calculation.getMSRLoudspeakerSignals( ...
    Input_Signal, setup, signal_info, system_info );

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