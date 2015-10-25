function [ setup ] = createSetup( settings )
%CREATESETUP Summary of this function goes here
%   Detailed explanation goes here


end

function setup = MultizoneSoundfieldSetup(res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi, phiL)
%%
c = 343;
k = (fmax/c)*2*pi;
M = ceil(k*R);

%%
quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, rq, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, rb, 'pw', 1.0, PWangb);
quiet.res  = res;
bright.res = quiet.res;
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');
bright = bright.setDesiredSoundfield(true, 'suppress_output');

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  disq, angq);
soundfield = soundfield.addSpatialZone(bright, disb, angb);

soundfield.BrightZ_Weight     = Wb;
soundfield.QuietZ_Weight      = Wq;
soundfield.UnattendedZ_Weight = Wu;

if N < 1
    N = 2*M + 1;
end
soundfield = soundfield.setN( N );

%%
setup = Speaker_Setup.loudspeaker_setup;
setup = setup.addMultizone_Soundfield(soundfield);
setup.Speaker_Arc_Angle = phiL;
setup.Angle_FirstSpeaker = phi;

if L < 1
    L = ceil( setup.Speaker_Arc_Angle/360 * 2*M + 1 );
end
setup.Loudspeaker_Count = L;

setup = setup.setRadius( Rl );

end