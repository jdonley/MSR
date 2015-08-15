function [ RIR ] = RIR_from_loudspeaker_setup_MCRoomSim( loudspeaker_setup )
%RIR_FROM_LOUDSPEAKER_SETUP Returns the RIR for each sample point in each
%zone for a given multizone loudspeaker setup using MCRoomSim
%   Returns an matrix which contains the RIR (column values) for each
%   sample point (rows).
%
%
%   Author: Jacob Donley, University of Wollongong, Australia
%   Email: Jacob.Donley089@uowmail.edu.au
%

%Set up room dimensions and characteristics
room = SetupRoom('Dim', [4.8, 3.3, 3]);
reproduction_center = [2.4 1.65 1.5];

%Set MCRoomSim advanced options here 
opts = MCRoomSimOptions;

%Add all sources (loudspeaker locations)
%s1 = AddSource('Location',[ 1 1 1 ]);
%s1 = [s1 AddSource('Location',[ 1 2 1 ])];
%s1 = [s1 AddSource('Location',[ 2 2 1 ])];
for spkr = 1:loudspeaker_setup.Loudspeaker_Count
   spkr_loc = loudspeaker_setup.Loudspeaker_Locations(spkr,:);
   [x, y] =  pol2cart(spkr_loc(1),spkr_loc(2));
   z = 0;
   src_loc = [x y z] + reproduction_center;
   src(spkr) = AddSource('Location', src_loc);
end

%Add all receviers (multizone sample point locations)
%r1 = AddReceiver('Location',[3 3 1]);
%r1 = [r1 AddReceiver('Location',[2 4 1])];
pt_good=1;
for pt = 1:numel(loudspeaker_setup.Bright_Samples)
    sample_loc = loudspeaker_setup.Bright_Samples_Locations(pt,:);
    if (~isnan(sample_loc))
        [x, y] =  pol2cart(sample_loc(1),sample_loc(2));
        z = 0;
        rec_loc = [x y z] + reproduction_center;
        rec(pt_good) = AddReceiver('Location', rec_loc);
        pt_good = pt_good+1;
    end
end
Bright_pts = pt_good;
pt_good=0;
for pt = 1:numel(loudspeaker_setup.Quiet_Samples)
    sample_loc = loudspeaker_setup.Quiet_Samples_Locations(pt,:);
    if (~isnan(sample_loc))
        [x, y] =  pol2cart(sample_loc(1),sample_loc(2));
        z = 0;
        rec_loc = [x y z] + reproduction_center;
        rec(Bright_pts + pt_good) = AddReceiver('Location', rec_loc);
        pt_good = pt_good+1;
    end
end

%Set all sources' settings
src = arrayfun(@(s) setfield(s,'Fs',16000),src);

%Set all receivers' settings
rec = arrayfun(@(s) setfield(s,'Fs',16000),rec);

%%
%Evalute the Room Impulse Responses
rir = RunMCRoomSim(src,rec,room,opts);

%%
rir_length = length(rir{1,1});
rir_temp = cell2mat(rir);
rir_temp = sum(rir_temp,2);
RIR = reshape(rir_temp,rir_length,length(r1));

end

