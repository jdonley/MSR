function ParametricMask( Signal_Length, LUT_resolution, Noise_Mask_dB, setup, mutlizone_setup )
%Clean_from_LUT_ZoneWeightedMask Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables
Drive = 'Z:\';

signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.L_noise_mask = Noise_Mask_dB; % dB
signal_info.method = 'ParametricMasker';
signal_info.weight = 1;
signal_info.input_filename = 'Masker';

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, LUT_resolution, Drive, 'new');


loudspeakers   = setup.Loudspeaker_Count;

% Generate Input Signal
level_mag = db2mag(signal_info.L_noise_mask);
%Input_Signal = v_addnoise( zeros(Signal_Length,1), signal_info.Fs, -Inf); % White noise
%Input_Signal = Perceptual_Tools.GreyNoise( Signal_Length/signal_info.Fs, signal_info.Fs, level_mag ); % Equal loudness shaped noise (ISO226)
Input_Signal = v_addnoise( zeros(Signal_Length,1), signal_info.Fs, -Inf, '', 5 ); % Speech shaped noise (SSN) (ITU-T P.50 Spectrum)

%% Shape noise spectrum to match quiet zone leakage spectrum
[spect,frqs] = Soundfield_Database.getQuietZoneSpectrumFromDB( mutlizone_setup, LUT_resolution, Drive, signal_info.Fs );
Input_Signal = Tools.shapeSpectrum( Input_Signal, spect, frqs, signal_info.Fs );


%% High-Pass Filter for speech reproduction
f_cutoff = signal_info.f_low;
% Design HIGH pass filter
[b,a] = butter(6, f_cutoff./(signal_info.Fs/2), 'high');
%Apply low pass filter to noise
Input_Signal_filt = filter(b,a,Input_Signal(:));
%Power normalise
Input_Signal = Broadband_Tools.power_norm( Input_Signal(:), Input_Signal_filt(:), signal_info.Fs, [f_cutoff, signal_info.Fs/2]);



%% Loudspeaker signals
Loudspeaker_Signals = repmat(Input_Signal(:),[1 loudspeakers]);


% Scale signals so they don't clip upon saving
if signal_info.L_noise_mask <= 0
    scaler = 1 / (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
elseif signal_info.L_noise_mask > 0 % For a positive masker we scale the signals to save up to a 40db noise masker
    scaler = 1 / (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
end

% Normalise Loudspeaker Signals
Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)))  ...
    *  scaler  ...
    *  level_mag;



%% Once we have the speaker signals we should save them for later use as .wav files
Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    [], [], [], [], ...
    signal_info.Fs );


end

