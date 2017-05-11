
noise_levels_vec = [...
    -Inf, ...   % Clean Speech is -Inf dB noise
    -10, ...    % Good Quality and Low Intelligibility is -10 dB noise
    ];


SYS = Current_Systems.IEEETransactions_System_C; % This is good for circular or line array and Quality/Privacy tuned noise

% Playback master gain
master_gain = -50; %dB (usually -50dB for recordings)








%%
% If a realworld recording is not specified in the system then abort
if ~any(strcmpi(strrep(SYS.signal_info.recording_type,'-',''),'realworld')), delete(gcp('nocreate')); return; end

%%
%Flip loudspeaker order (effectively flips entire setup) (false if not needed)
% (Only do this if loudspeakers have been calibrated in reverse order too)
MirrorSetup = SYS.system_info.MirrorSetup;

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
SYS.signal_info.L_noise_mask = -inf;
masker_signal_info = SYS.signal_info;

%%
paired = isfield(SYS.signal_info,'methods_list_paired') && SYS.signal_info.methods_list_paired;

ml_tmp = SYS.signal_info.methods_list;
for c = SYS.signal_info.methods_list_clean
    masker_list = SYS.signal_info.methods_list_masker;
    if paired
        masker_list = masker_list(c);
    end
    for m = masker_list
        m(m<1)=[];
        
        SYS.signal_info.method = [ml_tmp{c}];
        masker_signal_info.method = [ml_tmp{m}];
        
        % If we are performing an SPL analysis then there is no masker and
        % the masker method is used for the recording path as NoMask
        if isempty(masker_signal_info.method) && any(strcmp(SYS.analysis_info.Measures,'SPL'))
            masker_signal_info.method = SYS.signal_info.method;
        end
        
        subSYS = SYS;
        subSYS.Main_Setup(~(c==SYS.signal_info.methods_list_clean))=[];
        subSYS.Masker_Setup(~(m==SYS.signal_info.methods_list_masker))=[];
        
        if strcmpi( ...
                SYS.system_info.CurrentSpeakerArrayType, ...
                subSYS.Main_Setup.Speaker_Array_Type)
            for noise_mask = noise_levels_vec
                
                % masker_signal_info = [];
                % SYS.Masker_Setup = [];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                masker_signal_info.L_noise_mask = noise_mask; % dB
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %% Playback from setup, setup_info and system_info
                a = Hardware_Control.playMSR_andRec( subSYS.Main_Setup, SYS.Room_Setup, SYS.signal_info, SYS.system_info, subSYS.Masker_Setup, masker_signal_info, master_gain );
                
                fprintf('\nFinished noise mask level %d \n\n',noise_mask);
                
            end
        end
    end
end
SYS.signal_info.methods_list = ml_tmp;
