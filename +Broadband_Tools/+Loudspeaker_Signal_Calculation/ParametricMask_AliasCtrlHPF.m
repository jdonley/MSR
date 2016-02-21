function ParametricMask_AliasCtrlHPF( Signal_Length, LUT_resolution, Noise_Mask_dB, setup, mutlizone_setup )
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
signal_info.method = 'ParametricMaskerAliasCtrlHPF';
signal_info.weight = 1;
signal_info.input_filename = 'Masker';

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, LUT_resolution, Drive, 'new');


%% Generate Input Signal
level_mag = db2mag(signal_info.L_noise_mask);
Input_Signal = Perceptual_Tools.GreyNoise( Signal_Length/signal_info.Fs, signal_info.Fs, level_mag );

%% High-Pass Filter the noise with cutoff equivalent to multizone aliasing frequency
% This method formulates the cutoff frequency where aliasing begins to
% occur. It is based on the frequency where aliasing occurs for the
% multizone reproduction region but not necessarily the individual zones.
BZ = mutlizone_setup.Multizone_Soundfield.Bright_Zone;
QZ = mutlizone_setup.Multizone_Soundfield.Quiet_Zone;

R_ = max( [BZ.Radius_q + BZ.Origin_q.Distance; ...
    QZ.Radius_q + QZ.Origin_q.Distance;  ]);
phiL_rad = mutlizone_setup.Speaker_Arc_Angle / 180 * pi;

f_cutoff = signal_info.c * (mutlizone_setup.Loudspeaker_Count - 1) / (2 * R_ * phiL_rad);

% Try a higher cutoff for the parametric (one octave higher)
f_cutoff = f_cutoff * 2;

% Design HIGH pass filter
[b,a] = butter(6, f_cutoff./(signal_info.Fs/2), 'high');

%Apply low pass filter to noise
Input_Signal_filt = filter(b,a,Input_Signal(:));

%Adjust the power of the signal in the passband to match the non-filtered version
%(power normalisation (equalisation))
input_norm_sig = Input_Signal(:) ./ max(abs(Input_Signal(:)));
Input_Signal = Broadband_Tools.power_norm( input_norm_sig, Input_Signal_filt(:), signal_info.Fs, [f_cutoff, signal_info.Fs/2]);

%% Loudspeaker Signals
Loudspeaker_Signals = repmat(Input_Signal(:),[1 setup.Loudspeaker_Count]);

%% Scale signals so they don't clip upon saving
if signal_info.L_noise_mask <= 0
    scaler = 1 / (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
elseif signal_info.L_noise_mask > 0 % For a positive masker we scale the signals to save up to a 40db noise masker
    scaler = 1 / (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
end

%% Apply scaling and magnitude
Loudspeaker_Signals = Loudspeaker_Signals  ...
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

