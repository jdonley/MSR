clear;
clc;
clear;
clear classes;
close all;
tic;

%%
%%%% NEW LAYOUT %%%%
speech_layout = {'brightzone_pos_angle',        90, ...
    'quietzone_pos_angle',         -90, ...
    'brightzone_source_angle',     0};

array_type = 'circle';
RIRs = [];

for spkr = 1:2
delete(gcp);
    if strcmp(array_type, 'circle')
        [x,y] = pol2cart(-90/180*pi, 0.6);
        x_ = sqrt(1.5^2-y^2);
        th_c = atan2(y,-x_);
        th = th_c;
    elseif strcmp(array_type, 'line')
        x_=1.5;
        th_c = 180;
        th = atan2(-0.6,-1.5);
    end
    if spkr ==1
        %regular
        array = {'angleto_firstloudspeaker',      90, ...
            'numberof_loudspeakers',         24, ...
            'loudspeaker_model',             'Genelec 8010A', ...
            'loudspeaker_radius',            1.5, ...
            'loudspeaker_spacing',           [], ...
            'speaker_array_type',            array_type};
    elseif spkr == 2
        %parametric
        array = {'angleto_firstloudspeaker',     th/pi*180, ...
            'loudspeaker_radius',           x_, ...
            'numberof_loudspeakers',        1, ...
            'loudspeaker_model',            'Parametric', ...
            'loudspeaker_spacing',           0.01, ...
            'speaker_array_type',            'line'};
    end
    
    setup = Speaker_Setup.createSetup({...
        speech_layout{:}, ...
        array{:}, ...
        'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
        'reproduction_radius',          1.0, ...
        'bright_weight',                1.0, ...
        'quiet_weight',                 0, ...
        'unattended_weight',            0.05, ...
        'brightzone_source_dist',       x_, ...
        'brightzone_source_type',       'ps', ...
        'brightzone_radius',            0.3, ...
        'brightzone_pos_distance',      0.6, ...
        'quietzone_radius',             0.3, ...
        'quietzone_pos_distance',       0.6, ...
        'maximum_frequency',            8000, ...
        'angleof_loudspeakerarc',       180, ...
        'angleof_loudspeakerarrcentre',	180, ...
        'loudspeaker_object',           Parametric_Synthesis.parametric_soundfield });
    
    setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
    setup = setup.calc_Loudspeaker_Weights();
    setup = setup.reproduceSoundfield('DEBUG');
    
    %% RIR Generation for a particular setup...
    room = Room_Acoustics.Room;
    room.Room_Size = [10 10 10];
    %room.Room_Size = [4 9 3]; %35.G46e
    %room.Room_Size = [8 10 3]; %6.107
    %room.Room_Size = [9 14 3]; % Out to lunch
    
    room.Room_Dimensions = 3;
    room.Reproduction_Centre = room.Room_Size ./ 2;
    
    room.Wall_Absorb_Coeff = 1.0;
    %room.Wall_Absorb_Coeff = 0.3;
    
    reflect_coeff = sqrt(1-room.Wall_Absorb_Coeff);
    reverb_time = repmat(reflect_coeff,1,6); %Seconds
    Fs = 16000;
    n_samples = 8000;
    
    room.NoReceivers = 32;
        
    if spkr ==1
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
            setup, ...
            room, ...
            reverb_time, ...
            n_samples);
        
    RIRs = struct('Bright_RIRs', RIR_B, ...
        'Bright_Receiver_Positions', Rec_Bright_Pos, ...
        'Quiet_RIRs', RIR_Q, ...
        'Quiet_Receiver_Positions', Rec_Quiet_Pos);
    
    elseif spkr == 2
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos ] = Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
            setup, ...
            room, ...
            reverb_time, ...
            n_samples, ...
            RIRs);
        
    RIRs = struct('Bright_RIRs', RIR_B, ...
        'Bright_Receiver_Positions', Rec_Bright_Pos, ...
        'Quiet_RIRs', RIR_Q, ...
        'Quiet_Receiver_Positions', Rec_Quiet_Pos, ...
        'Matched_Receivers', {RIR_DB_fullpath});
    end
    
    
    
    %% Group delay
    %Aim the reproduction at the centre of the bright zone using a group delay
    %for all speakers (no delay implies the reproduction is aimed at the centre
    %of the loud speaker array)
    % C = 343; %Speed of sound m/s
    % max_delay = ceil(soundfield.Radius*2 / C * Fs);
    %
    % [Spkr_LocX, Spkr_LocY] = pol2cart( ...
    %                          setup.Loudspeaker_Locations(:,1), ...
    %                          setup.Loudspeaker_Locations(:,2));
    % Spkr_Loc = [Spkr_LocX, Spkr_LocY];
    %
    % Bright_Loc = [repmat(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X,size(Spkr_LocX)), ...
    %               repmat(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y,size(Spkr_LocY))];
    %
    % Delay_Dist = sum((Bright_Loc - Spkr_Loc) .^ 2, 2) .^ 0.5; %Euclidean distance
    %
    % Time_Delay_Samples = round( (max(Delay_Dist) - Delay_Dist) /C * Fs ) + 1;
    %
    % Time_Delay_TFs = zeros(setup.Loudspeaker_Count, max_delay);
    % for i = 1:setup.Loudspeaker_Count
    %     Time_Delay_TFs(i,Time_Delay_Samples(i))=1;
    % end
    
    
    %% Save the RIRs to a database for reuse
    RIR_DB_fullpath = Room_Acoustics.getRIRDatabasePath( ...
        setup, ...
        room, ...
        'Z:\', ...
        'new2');
    RIRDBpath = fileparts(RIR_DB_fullpath);
    
    if ~exist( RIRDBpath,'dir'); mkdir( RIRDBpath ); end
    save( RIR_DB_fullpath, ...
        'RIRs');
    %'Time_Delay_TFs');
    
    if spkr == 2
        DBPath = RIRs.Matched_Receivers;
        DB = load( DBPath );
        DB.RIRs.Matched_Receivers = RIR_DB_fullpath;
        RIRs = DB.RIRs;
        save( DBPath, 'RIRs');
    end
    
    %%
    %clear
    %
    hold on;
    scatter(rec_b(:,1),rec_b(:,2),'.g'); hold on
    scatter(rec_q(:,1),rec_q(:,2),'.y');
    scatter(Rec_Bright_Pos(:,1),Rec_Bright_Pos(:,2),'ob'); hold on;
    scatter(Rec_Quiet_Pos(:,1),Rec_Quiet_Pos(:,2),'or'); hold on;
    
    src = [];
    [src(:,1), src(:,2)] = pol2cart( setup.Loudspeaker_Locations(:,1), setup.Loudspeaker_Locations(:,2));
    src = [src zeros(size(src,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre, size(src,1),1);
    scatter(src(:,1),src(:,2),'ok');hold off;
    axis equal;
    %axis([0 room.Room_Size(1) 0 room.Room_Size(2)]);
    
end

%%
toc

    



