S = Parametric_Synthesis.parametric_soundfield;

S.res = 100;
S.field_size = [3 3];

S.Angle = 0;
S = S.set_fd(4000);

S.Origin = [0.0 0.9+tand(-S.Angle)*1.5 ];


S = S.createSoundfield;

%S.plotSoundfield((S.Soundfield),'complex'); hold on;

%contour3(abs(S.Soundfield),db2mag(-120:10:0),'-k');
%hold off;



setup.Soundfield_reproduced = S.Soundfield;