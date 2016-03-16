function playMSR( Main_Setup, signal_info, system_info, Masker_Setup, masker_signal_info )
%PLAYMSR Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    Masker_Setup = [];
    masker_signal_info=[];
end

%% Get signals
SpkrSignals = Hardware_Control.getMSR_Audio( Main_Setup, signal_info, system_info, Masker_Setup, masker_signal_info );


%% Playback
playrec('play', SpkrSignals, system_info.playbackChannels);

end

