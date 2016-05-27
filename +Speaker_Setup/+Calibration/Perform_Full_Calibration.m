clc;clear;
delete(gcp('nocreate'));
current_pool = parpool; %Start new pool

%% Settings
system_info.dev_model = 'ASIO Hammerfall DSP';
system_info.fs = 48000;
system_info.playbackChannels = ...
    [ 1  2  3  4  5  6  7  8 ...
    9 10 11 12 13 14 15 16 ...
    17 18 19 20 21 22 23 24];

system_info.recordChannels = ...
    [ 3 ];

system_info.Drive = 'Z:\';

fs = 48000;
f_low = 100;
f_high = 10000;
frequency_range = [f_low f_high];
Nspkrs = 24;
sweep_len = 10; %sec
sweep_endLengths = 1; %sec
EQ_filt_len = 0.5; %seconds
EQ_filt_reg = [60 -6]; %[passband_gain, stopband_gain] (dB)

Calibration_Path = [system_info.Drive '+Calibration_Data' filesep];
Filters_Dir = ['+Filters' filesep];
Recordings_Dir = ['+Recordings' filesep];

%% Generate Exponential Sine Sweep (ESS) and Inverse FFT
[y, invY] = Tools.synthSweep( ...
    sweep_len+sweep_endLengths, ...
    fs,...
    frequency_range(1),...
    frequency_range(2),...
    sweep_endLengths*fs);
y = [zeros(sweep_endLengths*fs,1); y(:)];

%% Duplicate Signals for Playback
ESS_len = length(y);
y_multi = reshape( ...
    [repmat( [y; zeros(ESS_len * Nspkrs,1)], Nspkrs - 1, 1); y], ...
            ESS_len * Nspkrs, ...
            Nspkrs);
        
%% Find and initialise hardware
if playrec('isInitialised')
    playrec('reset');
end

devs = playrec('getDevices');
for d = 1:length(devs)
    if strcmpi(devs(d).name,system_info.dev_model)
        dev = devs(d);
        break;
    end
end

playrec( 'init', system_info.fs, dev.deviceID, dev.deviceID );

%% Play and Record Each ESS
recID = playrec('playrec', ...
    y_multi, ...
    system_info.playbackChannels, ...
    -1, ... % Record the same number of samples as the playback
    system_info.recordChannels);

%% Wait for recording
fprintf('\n====== Playing and Recording Calibration Signals ======\n');
fprintf('\tCompletion: ');n=0;h=[];t=0;tic;

while ~playrec('isFinished',recID)
    pause(1); t=t+1;
    [n,h] = Tools.showTimeToCompletion(t/(size(y_multi,1)/system_info.fs), n, h);
end
fprintf('\n');

%% Retrieve recording
[y_multi_rec, recChanList] = playrec('getRec', recID );
if recChanList ~= system_info.recordChannels
    error('Recording completed with incorrect channel numbers');
end

%% Find Equalisation Filters
EQ = Speaker_Setup.Calibration.getCalibrationFilters( ...
    y_multi_rec, ...
    y, ...
    frequency_range, ...
    fs, ...
    EQ_filt_len, ...
    EQ_filt_reg );

%% Save Filters and Recordings
if ~exist([Calibration_Path Filters_Dir],'dir'); mkdir([Calibration_Path Filters_Dir]); end
save([Calibration_Path Filters_Dir 'EQ_Filters_' ...
    datestr(now,'yyyy-mm-dd_HH.MM') ... 
    '.mat'], 'EQ', 'fs');
if ~exist([Calibration_Path Recordings_Dir],'dir'); mkdir([Calibration_Path Recordings_Dir]); end
audiowrite([Calibration_Path Recordings_Dir 'Recording_' ...
    datestr(now,'yyyy-mm-dd_HH.MM') ... 
    '.wav'], y_multi_rec, fs);

%% Calibrate MSR Loudspeaker Signals
%Speaker_Setup.Calibration.Calibrate_MSR_Loudspeaker_Signals;





