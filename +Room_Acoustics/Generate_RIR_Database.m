clear;
clc;
%clear;
%close all;

%%
f_max = 8000;
quiet  = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'pw', 1.0, 15);
quiet.res  = 50;
bright.res = 50;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');    
bright = bright.setDesiredSoundfield(true);

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  0.60, 0);
soundfield = soundfield.addSpatialZone(bright, 0.60, 180);
%%
 soundfield.BrightZ_Weight     = 1.0;
 soundfield.QuietZ_Weight      = 1e4;
 soundfield.UnattendedZ_Weight = 0.05;

soundfield = soundfield.setN(80);
Radius = 1.0;
soundfield = soundfield.createSoundfield('DEBUG', Radius);


%%
 setup = Speaker_Setup.loudspeaker_setup;
 setup = setup.addMultizone_Soundfield(soundfield);
 %setup.Loudspeaker_Count = 2*ceil(f_max/343 * exp(1) * Radius / 2 ) + 1 ;%16;
 setup.Speaker_Arc_Angle = 360;
 setup.Angle_FirstSpeaker = 0;
 k = f_max/343*2*pi;
 M = ceil(k*Radius);
 setup.Loudspeaker_Count = ceil( setup.Speaker_Arc_Angle/360 * 2*M + 1 ) ;%16;
 setup = setup.setRadius(1.5);
 
 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

%% RIR Generation for a particular setup...

room_size = [5 6 4];
room_dimensions = 3;
reproduction_center = room_size ./ 2;
absorb_coeff = 0.3;
reflect_coeff = sqrt(1-absorb_coeff);
reverb_time = repmat(reflect_coeff,1,6); %Seconds
n_samples = 8000;
n_rec = 32;

[RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
    setup, ...
    room_size, ...
    room_dimensions, ...
    reproduction_center, ...
    reverb_time, ...
    n_samples, ...
    n_rec);

RIRs = struct('Bright_RIRs', RIR_B, ...
    'Bright_Receiver_Positions', Rec_Bright_Pos, ...
    'Quiet_RIRs', RIR_Q, ...
    'Quiet_Receiver_Positions', Rec_Quiet_Pos);


%% Save the RIRs to a database for reuse
RIR_Database_Path = '+Room_Acoustics\+RIR_Database\';

room = num2str(room_size(1));
room = [room 'x' num2str(room_size(2)) ];
if room_dimensions == 3
    room = [room 'x' num2str(room_size(3)) ];
end

room_cent = num2str(reproduction_center(1));
room_cent = [room_cent 'x' num2str(reproduction_center(2)) ];
if room_dimensions == 3
    room_cent = [room_cent 'x' num2str(reproduction_center(3)) ];
end

RIR_Name__Details = [num2str(setup.Loudspeaker_Count) 'Src_' num2str(n_rec) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(absorb_coeff) 'Ab'];
    
save([RIR_Database_Path 'RIRs__' RIR_Name__Details '.mat'],'RIRs');


%%
%clear
% 
scatter(rec_b(:,1),rec_b(:,2),'.g'); hold on
scatter(rec_q(:,1),rec_q(:,2),'.y');
scatter(Rec_Bright_Pos(:,1),Rec_Bright_Pos(:,2),'ob')
scatter(Rec_Quiet_Pos(:,1),Rec_Quiet_Pos(:,2),'or'); hold off;
axis equal;








