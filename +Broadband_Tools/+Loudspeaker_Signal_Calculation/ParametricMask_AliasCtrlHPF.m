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


% Generate Input Signal
level_mag = db2mag(signal_info.L_noise_mask);
Input_Signal = v_addnoise( zeros(Signal_Length,1), signal_info.Fs, -Inf); % White noise
%Input_Signal = Perceptual_Tools.GreyNoise( Signal_Length/signal_info.Fs, signal_info.Fs, level_mag ); % Equal loudness shaped noise (ISO226)
%Input_Signal = v_addnoise( zeros(Signal_Length,1), signal_info.Fs, -Inf, '', 5 ); % Speech shaped noise (SSN) (ITU-T P.50 Spectrum)

%% Shape the noise to the speech being used
[spect_sp,frqs_sp]=Tools.LTASS('M:\MSR\+Miscellaneous\+Speech_Files\');
% span=floor(length(frqs_sp(frqs_sp>=125 & frqs_sp<=8000))/7);
span = round(5/100 * length(spect_sp));

spect_sp = smooth( spect_sp, span);
Input_Signal = Tools.shapeSpectrum( Input_Signal, spect_sp, frqs_sp, signal_info.Fs );


%% Find aliasing frequency for use as cutoff
% This method formulates the cutoff frequency where aliasing begins to
% occur. It is based on the frequency where aliasing occurs for the
% reproduction region but not necessarily the individual zones.

BZ = mutlizone_setup.Multizone_Soundfield.Bright_Zone;
QZ = mutlizone_setup.Multizone_Soundfield.Quiet_Zone;

R_ = max( [BZ.Origin_q.Distance; ...
    	   QZ.Origin_q.Distance;  ]);
phiL_rad = mutlizone_setup.Speaker_Arc_Angle / 180 * pi;

f_cutoff = signal_info.c * (mutlizone_setup.Loudspeaker_Count - 1) / (2 * R_ * phiL_rad);


%% Shape noise spectrum to match quiet zone leakage spectrum
[spect,frqs] = Soundfield_Database.getQuietZoneSpectrumFromDB( mutlizone_setup, LUT_resolution, Drive, signal_info.Fs );
%spect = spect .^ 0.5; %square root to reduce emphasis of spectral adjustment
spect2 = spect;
span2 = round(5/100 * length(spect));

spect_low=spect(frqs<=f_cutoff);
db_drop = 50;
spect_low_new = linspace( ...
    db2mag(mag2db(spect_low(end))-db_drop), ...
    spect_low(end), ...
    length(spect_low) );
% spect_low_new = logspace( ...
%     log10(db2mag(mag2db(spect_low(end))-db_drop)), ...
%     log10(spect_low(end)), ...
%     length(spect_low) );
spect2(frqs<=f_cutoff) = spect_low_new;

spect = smooth( spect, span2); %Smooth arbitrary spectral adjustment
spect2 = smooth( spect2, span2); %Smooth arbitrary spectral adjustment

Input_Signal_toMatch = Tools.shapeSpectrum( Input_Signal, spect, frqs, signal_info.Fs );
Input_Signal = Tools.shapeSpectrum( Input_Signal, spect2, frqs, signal_info.Fs );


%% High-Pass Filter for speech reproduction
f_cutoff_sp = signal_info.f_low;
% Design HIGH pass filter
[b,a] = butter(6, f_cutoff_sp./(signal_info.Fs/2), 'high');
%Apply low pass filter to noise
Input_Signal_filt = filter(b,a,Input_Signal(:));
%Power normalise
Input_Signal = Broadband_Tools.power_norm( Input_Signal(:), Input_Signal_filt(:), signal_info.Fs, [f_cutoff_sp, signal_info.Fs/2]);

%% High-Pass Filter the noise with cutoff equivalent to multizone aliasing frequency
% Try a higher cutoff for the parametric (one octave higher)
 %f_cutoff = f_cutoff * 3;

% Design HIGH pass filter
[b,a] = butter(6, f_cutoff./(signal_info.Fs/2), 'high');

%Apply low pass filter to noise
Input_Signal_filt = filter(b,a,Input_Signal(:));

%Adjust the power of the signal in the passband to match the non-filtered version
%(power normalisation (equalisation))
input_norm_sig = Input_Signal_toMatch(:) ./ max(abs(Input_Signal_toMatch(:)));
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

