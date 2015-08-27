function [ RIR_Bright, RIR_Quiet, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = RIR_from_loudspeaker_setup_rir_generator( loudspeaker_setup, room_size, room_dimensions, reproduction_center, reverb_time, n_samples, n_rec )
%RIR_FROM_LOUDSPEAKER_SETUP Returns the RIR for each sample point in each
%zone for a given multizone loudspeaker setup using rir_generator
%   Returns an matrix which contains the RIR (column values) for each
%   sample point (rows).
%
%
%   Author: Jacob Donley, University of Wollongong, Australia
%   Email: Jacob.Donley089@uowmail.edu.au
%
tic;

%Set up room dimensions and characteristics
if nargin < 2
    room_size = [5, 4, 3];
end
if nargin < 3
    dim = 3; % Room dimension
else
   dim = room_dimensions;
end
if nargin < 4
    reproduction_center = room_size ./ 2;
end
if nargin < 5
    beta = 0.4; % Reverberation time (s)
else
    beta = reverb_time;
end
if nargin < 6
    n = 4096;   % Number of samples
else
    n = n_samples;
end
if nargin < 7
    n_rec = -1;   % Number of samples
end
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

%Set MCRoomSim advanced options here 
c = 340;    % Speed of sound (m/s)
Fs = 16000; % Sample frequency (samples/s) 

%Add all sources (loudspeaker locations)
src = [];
[src(:,1), src(:,2)] = pol2cart( loudspeaker_setup.Loudspeaker_Locations(:,1), loudspeaker_setup.Loudspeaker_Locations(:,2));
src = [src zeros(size(src,1),size(room_size,2)-2)] + repmat(reproduction_center, size(src,1),1);


%Add all receviers (multizone sample point locations)
%r1 = AddReceiver('Location',[3 3 1]);
%r1 = [r1 AddReceiver('Location',[2 4 1])];

mask_b = loudspeaker_setup.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask;
mask_q = loudspeaker_setup.Multizone_Soundfield.Quiet_Zone.Soundfield_d_mask;

rec_b = [];
angl = loudspeaker_setup.Bright_Samples_Locations(:,:,1);
dis = loudspeaker_setup.Bright_Samples_Locations(:,:,2);
[rec_b(:,1), rec_b(:,2)] = pol2cart( angl(:), dis(:));
rec_b = [rec_b zeros(size(rec_b,1),size(room_size,2)-2)] + repmat(reproduction_center, size(rec_b,1),1);
rec_b = rec_b(mask_b(:) ~= 0,:);

rec_q = [];
angl = loudspeaker_setup.Quiet_Samples_Locations(:,:,1);
dis = loudspeaker_setup.Quiet_Samples_Locations(:,:,2);
[rec_q(:,1), rec_q(:,2)] = pol2cart( angl(:), dis(:));
rec_q = [rec_q zeros(size(rec_q,1),size(room_size,2)-2)] + repmat(reproduction_center, size(rec_q,1),1);
rec_q = rec_q(mask_q(:) ~= 0,:);


Rec_Bright_Pos = rec_b;
Rec_Quiet_Pos = rec_q;
% Generate random number seed from current time
rng('shuffle');
%Use random samples
if n_rec ~= -1
    ind_rec_b_rnd = randi(size(Rec_Bright_Pos,1),1,n_rec);
    Rec_Bright_Pos = Rec_Bright_Pos(ind_rec_b_rnd,:);
    
    ind_rec_q_rnd = randi(size(Rec_Quiet_Pos,1),1,n_rec);
    Rec_Quiet_Pos = Rec_Quiet_Pos(ind_rec_q_rnd,:);
end   
% 
% scatter(rec_b(:,1),rec_b(:,2),'og'); hold on
% scatter(rec_q(:,1),rec_q(:,2),'og');
% scatter(Rec_Bright_Pos(:,1),Rec_Bright_Pos(:,2),'ob')
% scatter(Rec_Quiet_Pos(:,1),Rec_Quiet_Pos(:,2),'ob'); hold off;

%%
%Evalute the Room Impulse Responses
RIR_Bright = zeros([size(Rec_Bright_Pos,1), n, size(src,1)]);
RIR_Quiet  = zeros([size(Rec_Quiet_Pos,1), n, size(src,1)]);

current_pool = parpool; %Start new pool
parfor_progress( size(src,1) );
parfor s = 1:size(src,1)
    RIR_Bright(:,:,s) = rir_generator(c, Fs, Rec_Bright_Pos, src(s,:), room_size, beta, n , mtype, order, dim, orientation, hp_filter);
    
    RIR_Quiet(:,:,s) = rir_generator(c, Fs, Rec_Quiet_Pos, src(s,:), room_size, beta, n , mtype, order, dim, orientation, hp_filter);
        
	parfor_progress;
end
parfor_progress(0);

%%
tEnd = toc;
fprintf('\nRIR execution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script


end

