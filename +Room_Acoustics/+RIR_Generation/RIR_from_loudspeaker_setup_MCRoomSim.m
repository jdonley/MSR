function [ RIR_Bright, RIR_Quiet, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q, Fs_new ] = RIR_from_loudspeaker_setup_MCRoomSim( loudspeaker_setup, room, reverb_time, n_samples, rec_positions )
%RIR_FROM_LOUDSPEAKER_SETUP Returns the RIR for each sample point in each
%zone for a given multizone loudspeaker setup using MCRoomSim
%   Returns an matrix which contains the RIR (column values) for each
%   sample point (rows).
%
%
%   Author: Jacob Donley, University of Wollongong, Australia
%   Email: Jacob.Donley089@uowmail.edu.au
%
Fs = 16000;
if nargin < 2
    room.Room_Size = [10, 10, 10];
    room.Reproduction_Centre = room.Room_Size./2;
    room.NoReceivers = 8;
    room.Wall_Absorb_Coeff = 1.0;
end
if nargin < 3
    beta = 0.4; % Reverberation time (s)
else
    beta = reverb_time;
end
if nargin < 4
    n = Fs/2;   % Number of samples
else
    n = n_samples;
end
if nargin < 5
    rec_positions = [];
end
%Set up room dimensions and characteristics
MCroom = SetupRoom( ...
    'Dim', room.Room_Size, ...
    'Freq', [125, 250, 500, 1000, 2000, 4000], ...
    'Absorption', repmat(room.Wall_Absorb_Coeff,6,6) );

%Set MCRoomSim advanced options here 
MC_Fs = 48000;
MCopts = MCRoomSimOptions( ...
    'Verbose', true, ...
    'Fs', MC_Fs, ... %Minimum is 44.1kHz so we need to apply decimation later for lower sampling frequencies
    'Duration', n/Fs, ...
    'SoundSpeed', 343, ...
    'AutoCrop', false );

%% Add all sources (loudspeaker locations and orientations)
MCsrc = repmat( AddSource, [1 loudspeaker_setup.Loudspeaker_Count] );
%Locations
src = [];
[src(:,1), src(:,2)] = pol2cart( loudspeaker_setup.Loudspeaker_Locations(:,1), loudspeaker_setup.Loudspeaker_Locations(:,2));
src = [src zeros(size(src,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre, size(src,1),1);
%Orientations
src_dir_ang = loudspeaker_setup.Loudspeaker_Directions(:)/pi*180;
src_dir = [src_dir_ang, zeros(size(src_dir_ang,1),2)];


%% Compute individual directivity impulse responses 
MCnsamples = MC_Fs * 2;
[ImpRsp, DirLst] = Room_Acoustics.loudspeakerIR(loudspeaker_setup, 256, MCnsamples, MC_Fs);


%% Assign Source Values
%Assign differing values with loop
for spkr = 1:loudspeaker_setup.Loudspeaker_Count
    MCsrc(spkr).Location = src(spkr,:);
    MCsrc(spkr).Orientation = src_dir(spkr,:);
    MCsrc(spkr).Response = ImpRsp;
    MCsrc(spkr).Direction = DirLst(:,1:2);
end

%Set all similar sources' settings with arrayfun
MCsrc = arrayfun(@(s) setfield(s,'Type','impulse'),MCsrc);
MCsrc = arrayfun(@(s) setfield(s,'Fs',MC_Fs),MCsrc);

%% Add all receviers (multizone sample point locations)
if isempty(rec_positions)
    
    mask_b = loudspeaker_setup.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask;
    mask_q = loudspeaker_setup.Multizone_Soundfield.Quiet_Zone.Soundfield_d_mask;
    
    rec_b = [];
    angl = loudspeaker_setup.Bright_Samples_Locations(:,:,1);
    dis = loudspeaker_setup.Bright_Samples_Locations(:,:,2);
    [rec_b(:,1), rec_b(:,2)] = pol2cart( angl(:), dis(:));
    rec_b = [rec_b zeros(size(rec_b,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre, size(rec_b,1),1);
    rec_b = rec_b(mask_b(:) ~= 0,:);
    
    rec_q = [];
    angl = loudspeaker_setup.Quiet_Samples_Locations(:,:,1);
    dis = loudspeaker_setup.Quiet_Samples_Locations(:,:,2);
    [rec_q(:,1), rec_q(:,2)] = pol2cart( angl(:), dis(:));
    rec_q = [rec_q zeros(size(rec_q,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre, size(rec_q,1),1);
    rec_q = rec_q(mask_q(:) ~= 0,:);
    
    
    Rec_Bright_Pos = rec_b;
    Rec_Quiet_Pos = rec_q;
    % Generate random number seed from current time
    rng('shuffle');
    %Use random samples
    if room.NoReceivers ~= -1
        ind_rec_b_rnd = randi(size(Rec_Bright_Pos,1),1,room.NoReceivers);
        Rec_Bright_Pos = Rec_Bright_Pos(ind_rec_b_rnd,:);
        
        ind_rec_q_rnd = randi(size(Rec_Quiet_Pos,1),1,room.NoReceivers);
        Rec_Quiet_Pos = Rec_Quiet_Pos(ind_rec_q_rnd,:);
    end
    
else
    
    Rec_Bright_Pos = rec_positions.Bright_Receiver_Positions;
    Rec_Quiet_Pos = rec_positions.Quiet_Receiver_Positions;
    
end

MCrecb = repmat( AddReceiver, size(Rec_Bright_Pos,1), 1 );
MCrecq = repmat( AddReceiver, size(Rec_Quiet_Pos,1), 1 );
for receiver = 1:size(Rec_Bright_Pos,1)
    MCrecb(receiver).Location =  Rec_Bright_Pos(receiver,:);
    MCrecq(receiver).Location =  Rec_Quiet_Pos(receiver,:);
end

%Set all receivers' settings
MCrecb = arrayfun(@(s) setfield(s,'Fs',MC_Fs),MCrecb);
MCrecq = arrayfun(@(s) setfield(s,'Fs',MC_Fs),MCrecq);


%% Evalute the Room Impulse Responses
MCrirb = RunMCRoomSim(MCsrc,MCrecb,MCroom,MCopts);
MCrirq = RunMCRoomSim(MCsrc,MCrecq,MCroom,MCopts);

MC_n_samples = length(MCrirb{1,1});

RIR_Bright = permute( reshape( ...
    cell2mat(MCrirb(:)'), ...
    MC_n_samples, ...
    size(Rec_Bright_Pos,1), ...
    []), ...
    [2 1 3]);

RIR_Quiet = permute( reshape( ...
    cell2mat(MCrirq(:)'), ...
    MC_n_samples, ...
    size(Rec_Quiet_Pos,1), ...
    []), ...
    [2 1 3]);

% dec_rat = 2; % Decimation ratio
% n_ = ceil(size(RIR_Bright_,2)/dec_rat); % new length
% RIR_Bright = zeros([size(Rec_Bright_Pos,1), n_, size(src,1)]);
% RIR_Quiet  = zeros([size(Rec_Quiet_Pos,1), n_, size(src,1)]);
% for s = 1:size(RIR_Bright_,3)
%     for r = 1:size(RIR_Bright_,1)
%         RIR_Bright(r,:,s) = decimate(RIR_Bright_(r,:,s), dec_rat, 10);
%         RIR_Quiet(r,:,s)  = decimate(RIR_Quiet_(r,:,s), dec_rat, 10);
%     end
% end

Fs_new = MC_Fs;
end

