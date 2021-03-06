function [ RIR_Bright, RIR_Quiet, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = RIR_from_loudspeaker_setup_rir_generator( loudspeaker_setup, room, reverb_time, signal_info, rec_positions )
%RIR_FROM_LOUDSPEAKER_SETUP Returns the RIR for each sample point in each
%zone for a given multizone loudspeaker setup using rir_generator
%   Returns a matrix which contains the RIR (column values) for each
%   sample point (rows).
%
%
%   Author: Jacob Donley, University of Wollongong, Australia
%   Email: Jacob.Donley089@uowmail.edu.au
%
tic;

%Set up room dimensions and characteristics
if nargin < 2
    room.Room_Size = [5, 4, 3];
    room.Reproduction_Centre = room.Room_Size./2;
end
if nargin < 3
    beta = 0.4; % Reverberation time (s)
else
    beta = reverb_time;
end
if nargin < 4
    n = 4096;   % Number of samples
else
    n = signal_info.rir_duration * signal_info.Fs;
end
if nargin < 5
    rec_positions = [];
end

dim = room.Room_Dimensions;
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

c = signal_info.c;    % Speed of sound (m/s)
Fs = signal_info.Fs; % Sample frequency (samples/s) 

%Add all sources (loudspeaker locations)
src = [];
[src(:,1), src(:,2)] = pol2cart( loudspeaker_setup.Loudspeaker_Locations(:,1), loudspeaker_setup.Loudspeaker_Locations(:,2));
src = [src zeros(size(src,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre, size(src,1),1);


%Add all receviers (multizone sample point locations)
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

%%
%Evalute the Room Impulse Responses
RIR_Bright = zeros([size(Rec_Bright_Pos,1), n, size(src,1)]);
RIR_Quiet  = zeros([size(Rec_Quiet_Pos,1), n, size(src,1)]);

room_size = room.Room_Size;
current_pool = parpool; %Start new pool
fprintf('\n====== Building RIR Database ======\n');
parfor_progress( size(src,1) );
parfor s = 1:size(src,1)
    RIR_Bright(:,:,s) = rir_generator( ...
        c, Fs, ...
        Rec_Bright_Pos, ...
        src(s,:), ...
        room_size, ...
        beta, n , mtype, order, dim, orientation, hp_filter);
    
    RIR_Quiet(:,:,s) = rir_generator( ...
        c, Fs, ...
        Rec_Quiet_Pos, ...
        src(s,:), ...
        room_size, ...
        beta, n , mtype, order, dim, orientation, hp_filter);
        
	parfor_progress;
end
parfor_progress(0);

%%
tEnd = toc;
fprintf('\nRIR execution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script


end

