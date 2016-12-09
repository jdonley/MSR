
% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016
% Date: 20 April 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;
%close all

%fp = '+Miscellaneous\+Sine_Sweep\SineSweep_Log.WAV';
%fp = '+Miscellaneous\+TestAudio_Files\chirp_test.WAV';
%fp = '+Miscellaneous\+TestAudio_Files\sinusoids_150-8000Hz.WAV';
fp = '+Miscellaneous\+Speech_File_Test\Female1_SA2.WAV';

[ym, fs] = audioread(['M:\MSR\' fp]);

t_ = 0:1/fs:(length(ym)/fs - 1/fs);
ym = sin(2*pi*1000*t_);

% fs is sampling frequency of message
% fc is frequency of carrier signal
% fc_s is sampling frequency of carrier signal

ym = ym(:) ./ max(abs(ym(:)));
fc = 60000;

fc_s = 2*(fc+fs);
fc_s = ceil(fc_s/fs) * fs;

t = 0:1/fc_s:(length(ym)/fs - 1/fc_s);



P = db2mag(100); % Source pressure amplitude
m = 0.9; % Modulation index




%ym_up = interp1(ym,linspace(1,length(ym), length(ym)*ceil(fc*2/fs) ));
ym_up = upsample(ym, fc_s/fs )';

%ym_up = sin(2*pi*1000*t);% + .5*sin(2*pi*600*t) + 2*sin(2*pi*700*t);

ym_up_ = m * ym_up;

ym_up_hilb = hilbert(ym_up_);

% Square-Root Amplitude Modulation (SRAM)
E_1 = ( 1 + ym_up_ ).^(1/2);
y2_ = P * E_1 .* sin( fc * 2*pi * t);

% Single SideBand Amplitude Modulation using Hilbert Transform (SSB-AM)
E_ = (1 + ym_up_hilb);
y2 = P * real( E_ .* exp(1i * fc * 2*pi * t) );

E1 = real(E_1);
E  = real(E_);

%% Demodulation
% for d = 1:100
% distance=linspace(0.01,10,100);

%Beta = (gam+1)/2;
Beta = 1.2; % Coefficient of non-linearity
c = 343; %Speed of sound in air
rho = 1.225; % Density of medium
alpha = 0.7; %Absorption Coefficient of ultrasound carrier
A = 0.2; %Cross sectional area in metres^2
z = 1.5; % Axial distance in metres
tau = t - z/c;
%E = E of tau
%p1 = (Beta*P^2*A) / (16*pi*rho*c^4*z*alpha) * diff(diff(E1.^2)); %Berktay's Eq
%p2 = P * diff( diff( abs(real(E_)) ));
omega = 2 * pi * fc;% Angular frequency of carrier wave

p1 = (A * P) / (4*omega*pi*c*z) * diff(diff(  E1.*atan( Beta*omega*P*E1 / (4*alpha*rho*c^3) )  )); % Merklinger's Eq
p2 = (A * P) / (4*omega*pi*c*z) * diff(diff(  E .*atan( Beta*omega*P*E  / (4*alpha*rho*c^3) )  )); % Merklinger's Eq

% downsample
p1_down = downsample(p1, fc_s/fs);
p2_down = downsample(p2, fc_s/fs);
%normalise
p1_down = p1_down / max(abs(p1_down(:)));
p2_down = p2_down / max(abs(p2_down(:)));

% Nfft = 2^nextpow2(length(p1_down(:)));
% pow2(d) = max( mag2db(abs(fft(p2_down,Nfft))) );
% end
% figure(1)
% plot(distance,pow2);

%% Total Harmonic Distortion results
THD_perc1 = round(db2mag(thd(p1_down,fs,100))*100,4)
THD_perc2 = round(db2mag(thd(p2_down,fs,100))*100,4)

%% Sound result
soundsc(ym_up ,fc_s);
pause(3);

%soundsc(p2_down ,fs);
Tools.pesq_mex_vec(ym,p1_down,fs)
Tools.pesq_mex_vec(ym,p2_down,fs)

%% Plot result
% figure(10)
% plot(y2);       hold on;
% plot(ym_up); 
% plot(p_down);        hold off;

figure(20);
Nfft = 2^nextpow2(length(p1_down(:)));
Z1 = fft(p1_down,Nfft);
Z2 = fft(p2_down,Nfft);
Z_orig = fft(ym,Nfft);
f = linspace(0,fs/2,Nfft/2);
f = f(1:end-1);
plot(f, mag2db(abs( Z1(1:end/2-1) )) ,'r');  hold on;
plot(f, mag2db(abs( Z2(1:end/2-1) )) ,'b');  hold on;
plot(f, mag2db(abs( Z_orig(1:end/2-1) )) ,'k'); hold on;
hold off;
set(gca,'XScale','log');
%set(gca,'YScale','log');
grid on;
xlim([100 10000]);
ylim([-150 100]);
title('Original Message and Demodulated Message');


% figure(30);
% Nfft = 2^nextpow2(length(E_1(:)));
% Z_1 = fft(E_1,Nfft);
% Z_ = fft(E_,Nfft);
% f = linspace(0,fc_s/2,Nfft/2);
% f = f(1:end-1);
% plot(f, abs( Z_(1:end/2-1) ) ,'b'); hold on;
%  plot(f, abs( Z_1(1:end/2-1) ) ,'r'); hold off;
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% grid on;
% %xlim([fc, fc] + [-fs, +fs]);


% figure(40);
% Nfft = 2^nextpow2(length(y2_(:)));
% Z_1 = fft(y2_,Nfft);
% Z_ = fft(y2,Nfft);
% f = linspace(0,fc_s/2,Nfft/2);
% f = f(1:end-1);
% plot(f, mag2db(abs( Z_(1:end/2-1) )) ,'b'); hold on;
% plot(f, mag2db(abs( Z_1(1:end/2-1) )) ,'r');
% hold off;
% set(gca,'XScale','log');
% %set(gca,'YScale','log');
% grid on;
% xlim([fc, fc] + [-fs/2, +fs/2]);
% title('Modulated Signals');