clear;
clc;
clear;
clear classes;
%close all;
delete(gcp);
tic;

%%
loudspeaker_layout = {'numberof_loudspeakers',        24, ...
                      'loudspeaker_radius',           1.5, ...
                      'loudspeaker_model',            'Genelec 8010A', ...
                      'angleof_loudspeakerarrcentre', 180, ...
                      'loudspeaker_spacing',          0.01, ...
                      'speaker_array_type',           'line' };
                  
speech_layout = {'brightzone_pos_angle',        90, ...
                 'quietzone_pos_angle',         -90};
                  
% speech_layout = {'brightzone_pos_angle',        90, ...
%                  'quietzone_pos_angle',         -90, ...
%                  'brightzone_source_angle',     0};

Speech_Setup = Speaker_Setup.createSetup({ speech_layout{:}, loudspeaker_layout{:}});

Speech_Setup.Multizone_Soundfield = Speech_Setup.Multizone_Soundfield.createSoundfield('DEBUG');
Speech_Setup = Speech_Setup.calc_Loudspeaker_Weights();
Speech_Setup = Speech_Setup.reproduceSoundfield('DEBUG');

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

[RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
    Speech_Setup, ...
    room, ...
    reverb_time, ...
    n_samples);

RIRs = struct('Bright_RIRs', RIR_B, ...
    'Bright_Receiver_Positions', Rec_Bright_Pos, ...
    'Quiet_RIRs', RIR_Q, ...
    'Quiet_Receiver_Positions', Rec_Quiet_Pos);


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
RIRDBpath = fileparts(DB_fullpath);
    
if ~exist( RIRDBpath,'dir'); mkdir( RIRDBpath ); end
save( RIR_DB_fullpath, ...
        'RIRs');
        %'Time_Delay_TFs');


%%
%clear
% 
scatter(rec_b(:,1),rec_b(:,2),'.g'); hold on
scatter(rec_q(:,1),rec_q(:,2),'.y');
scatter(Rec_Bright_Pos(:,1),Rec_Bright_Pos(:,2),'ob'); hold on;
scatter(Rec_Quiet_Pos(:,1),Rec_Quiet_Pos(:,2),'or'); hold on;

src = [];
[src(:,1), src(:,2)] = pol2cart( Speech_Setup.Loudspeaker_Locations(:,1), Speech_Setup.Loudspeaker_Locations(:,2));
src = [src zeros(size(src,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre, size(src,1),1);
scatter(src(:,1),src(:,2),'ok');hold off;
axis equal;
%axis([0 room.Room_Size(1) 0 room.Room_Size(2)]);

%%
toc






