clc;
clear;

%%
SYS = Current_Systems.loadCurrentSRsystem;

%%
%Flip loudspeaker order (effectively flips entire setup) (false if not needed)
% (Only do this if loudspeakers have been calibrated in reverse order too)
MirrorSetup = false;

% Flip loudspeaker and mic order (effectively flips entire setup)
if MirrorSetup
    SYS.system_info.playbackChannels = flip(SYS.system_info.playbackChannels);
    SYS.system_info.recordChannels = flip(SYS.system_info.recordChannels);
end

% False or True to record reference signal
if SYS.signal_info.reference
    SYS.system_info.playbackChannels = SYS.signal_info.reference_channel;
end


%%
% Playback master gain
master_gain = 0; %dB

noise_levels_vec = SYS.signal_info.L_noise_mask;
noise_levels_vec(noise_levels_vec>0)=[];

SYS.signal_info.L_noise_mask = -inf;
masker_signal_info = SYS.signal_info;

%%
ml_tmp = SYS.signal_info.methods_list;
for c = SYS.signal_info.methods_list_clean
    for m = SYS.signal_info.methods_list_masker
        m(m<1)=[];
        
        SYS.signal_info.method = ml_tmp{c};
        masker_signal_info.method = ml_tmp{m};
        
        for noise_mask = noise_levels_vec
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            masker_signal_info.L_noise_mask = noise_mask; % dB
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %% Playback from setup, setup_info and system_info and record the result
            Hardware_Control.playMSR_andRec( SYS.Main_Setup, SYS.Room_Setup, SYS.signal_info, SYS.system_info, SYS.Masker_Setup, masker_signal_info, master_gain );
            
            fprintf('\nFinished noise mask level %d \n\n',noise_mask);
            
        end
    end
end
SYS.signal_info.methods_list = ml_tmp;
