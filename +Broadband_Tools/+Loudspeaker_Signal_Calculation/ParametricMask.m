function ParametricMask( Signal_Length, LUT_resolution, Noise_Mask_dB, weight, setup )
%Clean_from_LUT_ZoneWeightedMask Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables

signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.L_noise_mask = Noise_Mask_dB; % dB
signal_info.weight = weight;
signal_info.method = 'ParametricMask';

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, LUT_resolution, 'Z:\', 'new');


loudspeakers   = setup.Loudspeaker_Count;

% Generate Input Signal
level_mag = db2mag(Noise_Mask_dB);
Input_Signal = Perceptual_Tools.GreyNoise( Signal_Length/Fs, Fs, level_mag );
Loudspeaker_Signals = repmat(Input_Signal(:),[1 loudspeakers]);


% Scale signals so they don't clip upon saving
if Noise_Mask_dB <= 0
    scaler = 1 / (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
elseif Noise_Mask_dB > 0 % For a positive masker we scale the signals to save up to a 40db noise masker
    scaler = 1 / (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
end

% Normalise Loudspeaker Signals
Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)))  *  scaler;



%% Once we have the speaker signals we should save them for later use as .wav files
Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    loudspeakers, [], [], [], ...
    Fs );


end

