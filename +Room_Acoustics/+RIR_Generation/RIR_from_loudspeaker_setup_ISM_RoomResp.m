function [ RIR_Bright, RIR_Quiet ] = RIR_from_loudspeaker_setup_ISM_RoomResp( loudspeaker_setup )
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
room_size = [4.8, 3.3, 3];
reproduction_center = room_size ./ 2;
beta = 0.4; % Reverberation time (s)
n = 4096;   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 2;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

%Set MCRoomSim advanced options here 
c = 340;    % Speed of sound (m/s)
Fs = 16000; % Sample frequency (samples/s) 

%Add all sources (loudspeaker locations)
src = [];
[src(:,1), src(:,2)] = pol2cart( loudspeaker_setup.Loudspeaker_Locations(:,1), loudspeaker_setup.Loudspeaker_Locations(:,2));
src = [src zeros(size(src,1),1)] + repmat(reproduction_center, size(src,1),1);


%Add all receviers (multizone sample point locations)
%r1 = AddReceiver('Location',[3 3 1]);
%r1 = [r1 AddReceiver('Location',[2 4 1])];

rec_b = [];
[rec_b(:,1), rec_b(:,2)] = pol2cart( loudspeaker_setup.Bright_Samples_Locations(:,1), loudspeaker_setup.Bright_Samples_Locations(:,2));
rec_b = [rec_b zeros(size(rec_b,1),1)] + repmat(reproduction_center, size(rec_b,1),1);

rec_q = [];
[rec_q(:,1), rec_q(:,2)] = pol2cart( loudspeaker_setup.Quiet_Samples_Locations(:,1), loudspeaker_setup.Quiet_Samples_Locations(:,2));
rec_q = [rec_q zeros(size(rec_q,1),1)] + repmat(reproduction_center, size(rec_q,1),1);


%%
%Evalute the Room Impulse Responses
RIR_Bright = rir_generator(c, Fs, rec_b, src, room_size, beta, n ); %, mtype, order, dim, orientation, hp_filter);

RIR_Quiet = rir_generator(c, Fs, rec_q, src, room_size, beta, n ); %, mtype, order, dim, orientation, hp_filter);

%%
tEnd = toc;
fprintf('\nRIR execution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script


end

