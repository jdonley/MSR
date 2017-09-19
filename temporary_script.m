clear; clc;

c = 343;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [1.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = [1.5 1.5 1.5];              % Source position [x y z] (m)
L = [3 3 3];                % Room dimensions [x y z] (m)
n = 0.1*fs;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = 1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

betaA = (1 - [1 [1 1 1 1 1]*1]).^2;                 % Reverberation time (s)
hA = rir_generator(c, fs, r, s, L, betaA, n, mtype, order, dim, orientation, hp_filter);

beta1 = (1 - [0.00*1 [1 1 1 1 1]*1]).^2;                 % Reverberation time (s)
h1 = rir_generator(c, fs, r, s, L, beta1, n, mtype, order, dim, orientation, hp_filter);

% ri = [0.25 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
% si = [2.75 1.5 1.5];              % Source position [x y z] (m)
% hi = rir_generator(c, fs, ri, si, L, betaA, n, mtype, order, dim, orientation, hp_filter);

hf = h1(:)-hA(:);

rtxN = 24;
[yy,zz] = meshgrid(linspace(0,3,rtxN));
rtx = [zeros(rtxN*rtxN,1), yy(:), zz(:)];
% rtx = [0 1.5 1.5];
stx = s;              % Source position [x y z] (m)
htx = rir_generator(c, fs, rtx, stx, L, betaA, n, mtype, order, dim, orientation, hp_filter);
rrx = r;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
srx = rtx;              % Source position [x y z] (m)
for i = 1:size(rtx,1)
    hrx(i,:) = rir_generator(c, fs, rrx, srx(i,:), L, betaA, n, mtype, order, dim, orientation, hp_filter);
end

% htx = imag(hilbert(htx));

hc = Tools.fconv(htx.',hrx.');
hc = sum(hc(1:numel(hf),:),2);

% figure(1);
% % plot(hA); hold on;
% plot(hf); hold on
% % plot(hi); hold on;
% plot(hc); hold on;
% hold off

figure(2);
HF = fft(hf);
HC = fft(hc);

ff = linspace(0,fs/2,n/2+1)/1e3;ff(end)=[];

MagnitudeC = abs(HC);
MagnitudeC(end/2+1:end)=[];

PhaseDifference = mod(unwrap(angle(HF)) - unwrap(angle(HC)+pi),2*pi)/pi*180-180;
PhaseDifference(end/2+1:end)=[];

subplot(2,1,1);
plot(ff, mag2db(  MagnitudeC.*ff.'  ));
xlim([0.1 10]); %ylim([-60 0]);
grid on; grid minor; set(gca,'xscale','log');

subplot(2,1,2);
plot(ff,PhaseDifference,'-');
xlim([0.1 10]); ylim([-180 180]); 
grid on; grid minor; set(gca,'xscale','log');


