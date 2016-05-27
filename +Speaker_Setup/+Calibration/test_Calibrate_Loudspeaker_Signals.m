clc;clear;
delete(gcp);
current_pool = parpool; %Start new pool

%%
frequency_range = [50, 20000]; %Hz
regularisation = [60 -6]; %[passband_gain, stopband_gain] (dB)
filter_length = 0.5; %seconds

recording = audioread('Z:\+Recordings\Chamber\Rec2_4-3-16_11-30pm.wav');
[reference, fs] = audioread('Z:\+Recordings\Chamber\sinesweep_20Hz-20kHz_30sec_5secTail_5secHead.WAV');
%duration = length(audioread('Z:\+Speaker_Signals\CalibratedSweeps\Loudspeaker_1.wav'))/fs;


EQ = Speaker_Setup.Calibration.getCalibrationFilters( ...
    recording, ...
    reference, ...
    frequency_range, ...
    fs, ...
    filter_length, ...
    regularisation );

%% Example set
Nspkrs = 24;
fs_spkrSigs = 16000;
loudspeaker_signals = repmat(reference, [1, Nspkrs]);
loudspeaker_signals_upsampled = loudspeaker_signals;

% % Upsample to 48kHz
% up_factor = fs/fs_spkrSigs;
% [b,a] = butter( 9, 1/up_factor );
% loudspeaker_signals_upsampled = zeros(size(loudspeaker_signals,1)*up_factor, Nspkrs);
% for spkr = 1:Nspkrs
%  loudspeaker_signals_upsampled(:,spkr) = ...
%      filter( b, a, ...
%      interp( loudspeaker_signals(:,spkr), up_factor ) );
% end

% Calibrate
loudspeaker_signals_calibrated = Tools.fconv( ...
    loudspeaker_signals_upsampled, ...
    EQ );

% Save
loudspeaker_signals_calibrated = loudspeaker_signals_calibrated ./ max(abs(loudspeaker_signals_calibrated(:)));
for spkr = 1:Nspkrs
    %audiowrite(['Z:\+Speaker_Signals\CalibratedSweeps\Loudspeaker_Upsampled_' num2str(spkr) '.wav'], loudspeaker_signals_calibrated(:,spkr), fs);
    audiowrite(['Z:\+Speaker_Signals\CalibratedSweeps\Loudspeaker_' num2str(spkr) '.wav'], loudspeaker_signals_calibrated(:,spkr), fs);
end

save('Z:\+Speaker_Signals\CalibratedSweeps\EQ_Filters.mat', 'EQ', 'fs');


