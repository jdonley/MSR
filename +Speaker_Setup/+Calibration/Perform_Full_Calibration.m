function Perform_Full_Calibration(SYS, RecordTF)
% Author: Jacob Donley
% Email: jrd089@uowmail.edu.au

%% Load System
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end
if nargin < 2, RecordTF = false; end

% If a realworld recording is not specified in the system then abort
if ~any(strcmpi(strrep(SYS.signal_info.recording_type,'-',''),'realworld')), delete(gcp('nocreate')); return; end

% Playback master gain
master_gain = -20; %dB (usually -20dB for recordings)

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
% y_multi = reshape( ...
%     [repmat( [y; zeros(ESS_len * Nspkrs,1)], Nspkrs - 1, 1); y], ...
%     ESS_len * Nspkrs, ...
%     Nspkrs);
y_multi = repmat( y , 1, Nspkrs);

%%
for EQorTFloop = 1:(RecordTF+1)
    if EQorTFloop == 2
        y_multi = Speaker_Setup.Calibration.applyInverseFilters( ...
            y_multi, EQ ./ max(sqrt(sum(EQ.^2))) );
    end
    
    NCalibChans = numel(SYS.system_info.calibrationRecChannel);
    NRecChans = numel(SYS.system_info.recordChannels);
    
    y_multi_rec=[];
    for ch = 1:numel(SYS.system_info.playbackChannels)
        %% Play and Record Each ESS
        recID = playrec('playrec', ...
            y_multi(:,ch) .* db2mag(master_gain), ...
            SYS.system_info.playbackChannels(ch), ...
            -1, ... % Record the same number of samples as the playback
            [SYS.system_info.calibrationRecChannel, ...
            SYS.system_info.recordChannels]);
        
        % Wait for recording
        fprintf('\n====== Playing and Recording Calibration Signals ======\n');
        fprintf('\tCompletion: ');n=0;tic;
        
        while ~playrec('isFinished',recID)
            pause(1);
            n = Tools.showTimeToCompletion(toc/(size(y_multi,1)/SYS.system_info.fs), n);
        end
        fprintf('\n');
        
        % Retrieve recording
        [y_multi_rec(:,:,ch), recChanList] = playrec('getRec', recID );
        if recChanList ~= SYS.system_info.calibrationRecChannel
            error('Recording completed with incorrect channel numbers');
        end
        
    end
    
    %% Find Equalisation Filters and Get Transfer Functions
    [EQs,TFs] = Speaker_Setup.Calibration.getCalibrationFilters( ...
        y_multi_rec, y, SYS );
    
    if EQorTFloop == 1
        EQ = EQs(:,:,1:NCalibChans);
    elseif EQorTFloop == 2
        TF = TFs(:,:,(NCalibChans+1):(NCalibChans+NRecChans));
%         TF=TFs;
    end
end

%% Save Filters and Transfer Functions
fs = SYS.system_info.fs;
FD_dir = [SYS.system_info.Drive SYS.system_info.FilterData_dir];

if exist('EQ','var')
    if ~exist(FD_dir,'dir'); mkdir(FD_dir); end
    save([FD_dir 'EQ_Filters_' ...
        datestr(now,'yyyy-mm-dd_HH.MM') ...
        '.mat'], 'EQ', 'fs');
end

if exist('TF','var')
    if ~exist(FD_dir,'dir'); mkdir(FD_dir); end
    save([FD_dir 'Transfer_Functions_' ...
        datestr(now,'yyyy-mm-dd_HH.MM') ...
        '.mat'], 'TF', 'fs');
end

%% Save Recordings
y_multi_rec = reshape(y_multi_rec,size(y_multi_rec,1)*size(y_multi_rec,2),size(y_multi_rec,3));
CR_dir = [SYS.system_info.Drive SYS.system_info.CalibrationRec_dir];
if ~exist(CR_dir,'dir'); mkdir(CR_dir); end
audiowrite([CR_dir 'Recording_' ...
    datestr(now,'yyyy-mm-dd_HH.MM') ...
    '.wav'], y_multi_rec, SYS.system_info.fs, 'BitsPerSample', 64);

%%
fprintf('Calibration completed and filters saved.\n\n');




