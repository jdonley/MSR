function Recordings = playMSR_andRec( Main_Setup, Room_Setup, signal_info, system_info, Masker_Setup, masker_signal_info, gain_dB )
%PLAYMSR_ANDREC Summary of this function goes here

if nargin < 5
    Masker_Setup = [];
    masker_signal_info=[];
end

if nargin < 7
    gain_dB = 0;
end


%% Find and initialise hardware
if playrec('isInitialised')
    playrec('reset');
end

devs = playrec('getDevices');
dev = devs(strcmpi({devs.name}, system_info.dev_model));
if isempty(dev)
    wrnCol = [255,100,0]/255;
    cprintf(wrnCol, 'The hardware device model ''');
    cprintf(-wrnCol, [system_info.dev_model ' ']);fprintf('\b');
    cprintf(wrnCol, ''' was not found.\n');
    cprintf(wrnCol, 'Skipping play and record procedure.\n');
    return;
end

playrec( 'init', system_info.fs, dev.deviceID, dev.deviceID );


%% Get signals
[SpkrSignals, seg_details] = Hardware_Control.getMSR_Audio( Main_Setup, signal_info, system_info, Masker_Setup, masker_signal_info );


%% Playback and Record
recID = playrec('playrec', ...
    SpkrSignals .* db2mag(gain_dB), ...
    system_info.playbackChannels, ...
    -1, ... % Record the same number of samples as the playback
    system_info.recordChannels);

%% Wait for recording
fprintf('\n====== Playing and Recording Multizone Soundfield Reproduction ======\n');
fprintf('\tCompletion: ');n=0;h=[];t=0;tic;

while ~playrec('isFinished',recID)
    pause(1); t=t+1;
    [n,h] = Tools.showTimeToCompletion(t/(size(SpkrSignals,1)/system_info.fs), n, h);
end
fprintf('\n');

%% Retrieve recording
[Recordings, recChanList] = playrec('getRec', recID );

if recChanList ~= system_info.recordChannels
    error('Recording completed with incorrect channel numbers');
end


%% Split and Save recordings if not asked for as output argument
if nargout == 0
   Hardware_Control.splitSaveRecording( ...
       Recordings, ...
       seg_details, ...
       Main_Setup, ...
       Room_Setup, ...
       signal_info, ...
       masker_signal_info, ...
       system_info); 
end

end

