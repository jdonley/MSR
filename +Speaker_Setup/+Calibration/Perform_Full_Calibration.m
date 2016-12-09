clear;

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

% If a realworld recording is not specified in the system then abort
if ~any(strcmpi(SYS.signal_info.recording_type,'realworld')), delete(gcp('nocreate')); return; end

%% Find and initialise hardware
if playrec('isInitialised')
    playrec('reset');
end

devs = playrec('getDevices');
dev = devs(strcmpi({devs.name},SYS.system_info.dev_model));
if isempty(dev)
    wrnCol = [255,100,0]/255;
    cprintf(wrnCol, 'The hardware device model ''');
    cprintf(-wrnCol, [SYS.system_info.dev_model ' ']);fprintf('\b');
    cprintf(wrnCol, ''' was not found.\n');
    cprintf(wrnCol, 'Skipping calibration procedure.\n');
    delete(gcp('nocreate')); return;
end

playrec( 'init', SYS.system_info.fs, dev.deviceID, dev.deviceID );


%% Generate Exponential Sine Sweep (ESS) and Inverse FFT
[y, invY] = Tools.synthSweep( ...
    SYS.system_info.Sweep_Length+SYS.system_info.Sweep_EndBuffers, ...
    SYS.system_info.fs,...
    SYS.system_info.f_low,...
    SYS.system_info.f_high,...
    SYS.system_info.Sweep_EndBuffers*SYS.system_info.fs, ...
    -3 );
y = [zeros(SYS.system_info.Sweep_EndBuffers*SYS.system_info.fs,1); y(:)];

%% Duplicate Signals for Playback
ESS_len = length(y);
Nspkrs = numel(SYS.system_info.playbackChannels);
y_multi = reshape( ...
    [repmat( [y; zeros(ESS_len * Nspkrs,1)], Nspkrs - 1, 1); y], ...
            ESS_len * Nspkrs, ...
            Nspkrs);
        
%% Play and Record Each ESS
recID = playrec('playrec', ...
    y_multi, ...
    SYS.system_info.playbackChannels, ...
    -1, ... % Record the same number of samples as the playback
    SYS.system_info.calibrationRecChannel);

%% Wait for recording
fprintf('\n====== Playing and Recording Calibration Signals ======\n');
fprintf('\tCompletion: ');n=0;tic;

while ~playrec('isFinished',recID)
    pause(1);
    n = Tools.showTimeToCompletion(toc/(size(y_multi,1)/SYS.system_info.fs), n);
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
    [SYS.system_info.f_low, ...
    SYS.system_info.f_high], ...
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

%%
fprintf('Calibration completed and filters saved.\n\n');




