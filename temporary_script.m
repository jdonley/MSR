c = 343;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [1.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = [1.5 1.5 1.5];              % Source position [x y z] (m)
L = [3 3 3];                % Room dimensions [x y z] (m)
n = 0.05*fs;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = 2;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

betaA = (1 - [1 [1 1 1 1 1]*1]).^2;                 % Reverberation time (s)
hA = rir_generator(c, fs, r, s, L, betaA, n, mtype, order, dim, orientation, hp_filter);

beta1 = (1 - [0.00*1 [1 1 1 1 1]*1]).^2;                 % Reverberation time (s)
h1 = rir_generator(c, fs, r, s, L, beta1, n, mtype, order, dim, orientation, hp_filter);

ri = [0.25 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
si = [2.75 1.5 1.5];              % Source position [x y z] (m)
hi = rir_generator(c, fs, ri, si, L, betaA, n, mtype, order, dim, orientation, hp_filter);

hf = h1(:)-hA(:);

rtx = [0.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
stx = [1.5 1.5 1.5];              % Source position [x y z] (m)
htx = rir_generator(c, fs, rtx, stx, L, betaA, n, mtype, order, dim, orientation, hp_filter);
rrx = [1.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
srx = [0.0 1.5 1.5];              % Source position [x y z] (m)
hrx = rir_generator(c, fs, rrx, srx, L, betaA, n, mtype, order, dim, orientation, hp_filter);

hc = Tools.fconv(htx(:),hrx(:));

figure(1);
% plot(hA); hold on;
plot(hf); hold on
% plot(hi); hold on;
plot(hc); hold on;
hold off

figure(2);
HF = fft(hf);
HC = fft(hc);

plot(unwrap(angle(HF)-angle(HC)))