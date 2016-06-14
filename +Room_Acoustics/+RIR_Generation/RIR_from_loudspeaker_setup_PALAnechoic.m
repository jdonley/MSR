function [ RIR_Bright, RIR_Quiet, Rec_Bright_Pos, Rec_Quiet_Pos ] = RIR_from_loudspeaker_setup_PALAnechoic( loudspeaker_setup, room, signal_info, rec_positions, norm_positions )
%RIR_FROM_LOUDSPEAKER_SETUP Returns the RIR for each sample point in each
%zone for a given multizone loudspeaker setup using MCRoomSim
%   Returns a matrix which contains the RIR (column values) for each
%   sample point (rows).
%
%
%   Author: Jacob Donley, University of Wollongong, Australia
%   Email: Jacob.Donley089@uowmail.edu.au
%
n = signal_info.rir_duration * signal_info.Fs;

if nargin < 4
    rec_positions = [];
end
if nargin < 5
    norm_positions = [];
end


%% Compute individual directivity impulse responses
ImpRsp = Room_Acoustics.loudspeakerIR(loudspeaker_setup, room, 0, n, signal_info.Fs, rec_positions, norm_positions);

Rec_Bright_Pos = rec_positions.Bright_Receiver_Positions;
Rec_Quiet_Pos  = rec_positions.Quiet_Receiver_Positions;

B_Nrec = size(Rec_Bright_Pos,1);
Q_Nrec = size(Rec_Quiet_Pos,1);

RIR_Bright = ImpRsp(1:B_Nrec,:);
RIR_Quiet  = ImpRsp((B_Nrec+1):(B_Nrec + Q_Nrec),:);
end

