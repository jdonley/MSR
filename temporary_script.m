clear; clc;

c = 343;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
L = [3 3 3];                % Room dimensions [x y z] (m)
n = 0.2*fs;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = 1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter
rng shuffle;

betaA = (1 - [1 [1 1 1 1 1]*1]).^2;                 % Reverberation time (s)
beta1 = (1 - [0.00*1 [1 1 1 1 1]*1]).^2;                 % Reverberation time (s)

rtxN = 24;
[yy,zz] = meshgrid(linspace(0,3,rtxN));
rtx = [zeros(rtxN*rtxN,1), yy(:), zz(:)];
srx = rtx;

MM=[];PP=[];
tic;
ss=0;
while true %for ss = 1:10
    ss = ss+1;
% r = [1.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
r = rand(1,3)*3;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
% s = [1.5 1.5 1.5];              % Source position [x y z] (m)
s = rand(1,3)*3;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)

hA = rir_generator(c, fs, r, s, L, betaA, n, mtype, order, dim, orientation, hp_filter);
h1 = rir_generator(c, fs, r, s, L, beta1, n, mtype, order, dim, orientation, hp_filter);
hf = h1(:)-hA(:);

stx = s;              % Source position [x y z] (m)
htx = rir_generator(c, fs, rtx, stx, L, betaA, n, mtype, order, dim, orientation, hp_filter);
rrx = r;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
for i = 1:size(rtx,1)
    hrx(i,:) = rir_generator(c, fs, rrx, srx(i,:), L, betaA, n, mtype, order, dim, orientation, hp_filter);
end

% htx = imag(hilbert(htx));

hc = Tools.fconv(htx.',hrx.');
hc = sum(hc(1:numel(hf),:),2);
hc = Broadband_Tools.power_norm(hf,hc,fs,[250 1000]);

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

MM(:,ss) = MagnitudeC; %.*ff.';
PP(:,ss) = PhaseDifference;
% MM = mean([MM , MagnitudeC.*ff.' ],2);
% PP = mean([PP , PhaseDifference  ],2);


subplot(2,1,1);
plot(ff, mag2db(  MM  ),':k'); hold on;
plot(ff, mag2db(  mean(MM,2)  ),'-b','linew',1.5); hold off;
xlim([0.1 10]); %ylim([-60 0]);
grid on; grid minor; set(gca,'xscale','log');
xlabel('Frequency (kHz)');ylabel('Magnitude (dB)');

subplot(2,1,2);
plot(ff, PP ,':k'); hold on;
plot(ff, mean(PP,2) ,'-','co',.8*[1 0 0],'linew',1.5); hold off;
xlim([0.1 10]); ylim([-180 180]); yticks([-180:45:180]);
grid on; grid minor; set(gca,'xscale','log');
xlabel('Frequency (kHz)');ylabel('Phase (\circ)');

drawnow;
fprintf('%d\n',ss);
end
toc

