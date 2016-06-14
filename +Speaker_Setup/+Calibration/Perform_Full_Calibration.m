clc;clear;
delete(gcp('nocreate'));
current_pool = parpool; %Start new pool

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;


%% Generate Exponential Sine Sweep (ESS) and Inverse FFT
[y, invY] = Tools.synthSweep( ...
    SYS.system_info.Sweep_Length+SYS.system_info.Sweep_EndBuffers, ...
    SYS.system_info.fs,...
    frequency_range(1),...
    frequency_range(2),...
    SYS.system_info.Sweep_EndBuffers*SYS.system_info.fs);
y = [zeros(SYS.system_info.Sweep_EndBuffers*SYS.system_info.fs,1); y(:)];

%% Duplicate Signals for Playback
ESS_len = length(y);
Nspkrs = numel(SYS.system_info.playbackChannels);
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
    if strcmpi(devs(d).name,SYS.system_info.dev_model)
        dev = devs(d);
        break;
    end
end

playrec( 'init', SYS.system_info.fs, dev.deviceID, dev.deviceID );

%% Play and Record Each ESS
recID = playrec('playrec', ...
    y_multi, ...
    SYS.system_info.playbackChannels, ...
    -1, ... % Record the same number of samples as the playback
    SYS.system_info.calibrationRecChannel);

%% Wait for recording
fprintf('\n====== Playing and Recording Calibration Signals ======\n');
fprintf('\tCompletion: ');n=0;t=0;tic;

while ~playrec('isFinished',recID)
    pause(1); t=t+1;
    n = Tools.showTimeToCompletion(t/(size(y_multi,1)/SYS.system_info.fs), n);
end
fprintf('\n');

%% Retrieve recording
[y_multi_rec, recChanList] = playrec('getRec', recID );
if recChanList ~= SYS.system_info.calibrationRecChannel
    error('Recording completed with incorrect channel numbers');
end

%% Find Equalisation Filters
EQ = Speaker_Setup.Calibration.getCalibrationFilters( ...
    y_multi_rec, ...
    y, ...
    frequency_range, ...
    SYS.system_info.fs, ...
    SYS.system_info.Calibration_FiltLen, ...
    SYS.system_info.Calibration_FiltReg );

%% Save Filters and Recordings
fs = SYS.system_info.fs;

FD_dir = [SYS.system_info.Drive SYS.system_info.FilterData_dir];
if ~exist(FD_dir,'dir'); mkdir(FD_dir); end
save([FD_dir 'EQ_Filters_' ...
    datestr(now,'yyyy-mm-dd_HH.MM') ... 
    '.mat'], 'EQ', 'fs');

CR_dir = [SYS.system_info.Drive SYS.system_info.CalibrationRec_dir];
if ~exist(CR_dir,'dir'); mkdir(CR_dir); end
audiowrite([CR_dir 'Recording_' ...
    datestr(now,'yyyy-mm-dd_HH.MM') ... 
    '.wav'], y_multi_rec, SYS.system_info.fs);






