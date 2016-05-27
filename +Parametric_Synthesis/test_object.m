S = Parametric_Synthesis.parametric_soundfield;

S.res = 100;
S.field_size = [3 3];

S.Angle = 0;
S = S.set_fd(8000);

S.Origin = [1.5 1.5+tand(-S.Angle)*1.5 ];
S = S.set_f1( 36000 );
S = S.createEmptySoundfield();            
S = S.reproduceSoundfield('convolutional');

%S.plotSoundfield((S.Soundfield),'complex'); hold on;

%contour3(abs(S.Soundfield),db2mag(-120:10:0),'-k');
% hold off;


setup.Soundfield_reproduced = S.Soundfield;


%% Directivity plot
th = pi/2:-pi/180:-pi/2+pi/180;
r = repmat(1.5,size(th));
[x,y]=pol2cart(th,r*S.res);
O = length(S.Soundfield)/2;
response = abs(S.Soundfield((ceil(x+O)-1)*O*2 + ceil(y+O)));

figure(2); hold on
plot(th/pi*180,mag2db(response/max(response)))
xlim([-40,40]);ylim([-30,1]);grid on; grid minor; hold off;