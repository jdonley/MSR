function result = MSR(res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi, phiL)
fprintf('\n\nRunning MSR...\n\n');
%%
c = 343;
k = (f/c)*2*pi;
M = ceil(k*R);
kmax = (fmax/c)*2*pi;
Mmax = ceil(kmax*R);

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
soundfield.Radius = R;

soundfield.BrightZ_Weight     = Wb;
soundfield.QuietZ_Weight      = Wq;
soundfield.UnattendedZ_Weight = Wu;

if N < 1
    N = 2*M + 1;
end
soundfield = soundfield.setN( floor(N) );
soundfield = soundfield.createSoundfield('DEBUG', R);

%%
setup = Speaker_Setup.loudspeaker_setup;
setup = setup.addMultizone_Soundfield(soundfield);
setup.Speaker_Arc_Angle = phiL;
setup.Angle_FirstSpeaker = phi;

if L < 1
    L = ceil( setup.Speaker_Arc_Angle/360 * 2*Mmax + 1 );
end
setup.Loudspeaker_Count = floor(L);

setup = setup.setRadius( Rl );

 setup = setup.calc_Loudspeaker_Weights();
 setup = setup.reproduceSoundfield('DEBUG');

 
 %%
 figH = figure(...
      'name'	, 'Multizone Soundfield', ...
      'visible'	, 'on');
 realistic = false;
 setup.plotSoundfield( setup.Soundfield_reproduced, 'default', realistic);
 
 
 result = webfigure(figH);
 close(figH);
 
 end