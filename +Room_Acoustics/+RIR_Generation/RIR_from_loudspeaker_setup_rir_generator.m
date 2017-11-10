function [ RIR_Bright, RIR_Quiet, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = RIR_from_loudspeaker_setup_rir_generator( loudspeaker_setup, room, reverb_time, signal_info, rec_positions )
%RIR_FROM_LOUDSPEAKER_SETUP Returns the RIR for each sample point in each
%zone for a given multizone loudspeaker setup using rir_generator
%   Returns a matrix which contains the RIR (column values) for each
%   sample point (rows).
%
%

% Computes the RIR for each sample point in each zone for a given multizone loudspeaker setup using rir_generator
%
% Syntax:   [ RIR_Bright, RIR_Quiet, ...
%             Rec_Bright_Pos, Rec_Quiet_Pos, ...
%             rec_b, rec_q ] = ...
%                   RIR_from_loudspeaker_setup_rir_generator( ...
%                       loudspeaker_setup, room, reverb_time, ...
%                       signal_info, rec_positions )
%
% Inputs:
% 	loudspeaker_setup - The Speaker_Setup.loudspeaker_setup object
% 	room - The Room_Acoustics.Room object
% 	reverb_time - The reverberation time as described in rir_generator
% 	signal_info - The SR system objects signal_info structure
% 	rec_positions - (Optional) Structure of specific recevier positions
%                   (see Room_Acoustics.Generate_RIR_Database for example)
%
% Outputs:
% 	RIR_Bright - RIRs in the bright zone
% 	RIR_Quiet - RIRs in the quiet zone
% 	Rec_Bright_Pos - Position of the RIRs in the bright zone
% 	Rec_Quiet_Pos - Position of the RIRs in the quiet zone
% 	rec_b - All possible RIR bright zone positions
% 	rec_q - All possible RIR quiet zone positions
%
% Example:
% 	%See Room_Acoustics.Generate_RIR_Database
%
% See also: Generate_RIR_Database

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015-2017
% Date: 25 August 2017
% Version: 0.2 (25 August 2017)
% Version: 0.1 (15 August 2015)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


startTime = tic; %Start timing this function

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
if ~isempty(room.ReceiverPositions)
    rec_positions = struct('Bright_Receiver_Positions', ...
        room.ReceiverPositions(:,:,1), ...
        'Quiet_Receiver_Positions', ...
        nan(1,3));
    if size(room.ReceiverPositions,3) > 1
        rec_positions.Quiet_Receiver_Positions = ...
            room.ReceiverPositions(:,:,2);
    end
end

dim = room.Room_Dimensions;
mtype = room.ReceiverDirectivityPattern;    % Type of microphone
order = room.Reflection_Order;              % -1 equals maximum reflection order!
orientation = room.ReceiverOrientations;	% Microphone orientation [az el] (rad)
hp_filter = 0;                              % Enable high-pass filter

c = signal_info.c;    % Speed of sound (m/s)
Fs = signal_info.Fs; % Sample frequency (samples/s)

%Add all sources (loudspeaker locations)
src = loudspeaker_setup.Loudspeaker_Locations;
src = ...
    [src(:,1), ... % [azimuth]
    zeros(size(src,1),size(room.Room_Size,2)-size(src,2)), ...
    src(:,2:end)]; % [radius] or [elevation, radius]

[src(:,1), src(:,2), src(:,3)] = ...
    sph2cart( src(:,1), src(:,2), src(:,3));
src = src + repmat(room.Reproduction_Centre([2 1 3]), size(src,1),1);


%Add all receviers (multizone sample point locations)
% TODO: Edit the following so random receivers are also found in the 3rd dimension
if isempty(rec_positions) || ~isempty(room.ReceiverPositions)
    
    mask_b = loudspeaker_setup.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask;
    mask_q = loudspeaker_setup.Multizone_Soundfield.Quiet_Zone.Soundfield_d_mask;
    
    rec_b = [];
    angl = loudspeaker_setup.Bright_Samples_Locations(:,:,1);
    dis = loudspeaker_setup.Bright_Samples_Locations(:,:,2);
    [rec_b(:,1), rec_b(:,2)] = pol2cart( angl(:), dis(:));
    rec_b = [rec_b zeros(size(rec_b,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre([2 1 3]), size(rec_b,1),1);
    rec_b = rec_b(mask_b(:) ~= 0,:);
    
    rec_q = [];
    angl = loudspeaker_setup.Quiet_Samples_Locations(:,:,1);
    dis = loudspeaker_setup.Quiet_Samples_Locations(:,:,2);
    [rec_q(:,1), rec_q(:,2)] = pol2cart( angl(:), dis(:));
    rec_q = [rec_q zeros(size(rec_q,1),size(room.Room_Size,2)-2)] + repmat(room.Reproduction_Centre([2 1 3]), size(rec_q,1),1);
    rec_q = rec_q(mask_q(:) ~= 0,:);
end

% TODO: Edit the following so random receivers are also found in the 3rd dimension
if isempty(rec_positions)
    RmBounds = room.Room_Size([2 1 3]);
    % Only use positions that are within the room
    validPosInds = (rec_b(:,1) >= 0) & (rec_b(:,1) <= RmBounds(1)) & ...
        (rec_b(:,2) >= 0) & (rec_b(:,2) <= RmBounds(2)) & ...
        (rec_b(:,3) >= 0) & (rec_b(:,3) <= RmBounds(3));
    Rec_Bright_Pos = rec_b(validPosInds,:);
    Rec_Quiet_Pos  = rec_q(validPosInds,:);
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

if all(isnan(Rec_Bright_Pos(:)))
    Rec_Bright_Pos = nan*Rec_Quiet_Pos;
elseif all(isnan(Rec_Quiet_Pos(:)))
    Rec_Quiet_Pos = nan*Rec_Bright_Pos;
end

%% For boundary cancellation methods we will also synthesis an anechoic response

if any(contains(lower(signal_info.methods_list), 'boundarycancel' )) ... %if boundary cancel method exists in system
        && strcmpi(room.SystemType,'transmit') ...                       %and the room is set up to produce audio
        && loudspeaker_setup.Loudspeaker_Count == 1                      %and there is only one source of audio
    generateAnechoicRIRs = true;
else
    generateAnechoicRIRs = false;
end


%%
%Evalute the Room Impulse Responses
Nsrc = size(src,1);
Nrec = size(Rec_Bright_Pos,1);
Ntot = Nsrc*Nrec;
room_size = room.Room_Size([2 1 3]); % Needs to be [x,y,z]

RIR_Bright = zeros([Nrec n Nsrc (1+generateAnechoicRIRs)]); % initialise
RIR_Quiet = zeros([Nrec n Nsrc (1+generateAnechoicRIRs)]); % initialise

for anecho = 1:(1+generateAnechoicRIRs)
    if generateAnechoicRIRs && anecho == 2
        beta = beta*0;
    end
    
    RIR_Bright_ = zeros([Ntot, n]);
    RIR_Quiet_  = zeros([Ntot, n]);
    
    current_pool = gcp; %Start new pool
    fprintf('\n====== Building RIR Database ======\n');
    fprintf('\t Completion: '); startTime = tic;
    Tools.showTimeToCompletion;
    percCompl = parfor_progress( Ntot );
    parfor rs = 1:Ntot
        [r,s] = ind2sub([Nrec Nsrc],rs);
        RIR_Bright_(rs,:) = rir_generator( ...
            c, Fs, ...
            Rec_Bright_Pos(r,:), ...
            src(s,:), ...
            room_size, ...
            beta, n , mtype, order, dim, orientation, hp_filter);
        
        RIR_Quiet_(rs,:) = rir_generator( ...
            c, Fs, ...
            Rec_Quiet_Pos(r,:), ...
            src(s,:), ...
            room_size, ...
            beta, n , mtype, order, dim, orientation, hp_filter);
        
        %%%
        percCompl = parfor_progress;
        Tools.showTimeToCompletion( percCompl/100, [], [], startTime );
        %%%
    end
    
    % Reshape due to 1D parfor
    RIR_Bright(:,:,:,anecho) = permute( reshape(RIR_Bright_,[Nrec Nsrc n]), [1 3 2]);
    RIR_Quiet(:,:,:,anecho) = permute( reshape(RIR_Quiet_ ,[Nrec Nsrc n]), [1 3 2]);
    
    percCompl=parfor_progress(0);
    Tools.showTimeToCompletion( percCompl/100, [], [], startTime );
end
%%
tEnd = toc;
fprintf('\nRIR execution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script


end

