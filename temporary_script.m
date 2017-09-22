F = [linspace(0,.475,50) linspace(.525,1,50)];
H = [zeros(1,50) exp(-1j*pi*13*F(51:100))];
f = fdesign.arbmagnphase('Nb,Na,F,H',12,10,F,H);
W = [ones(1,50) 100*ones(1,50)];
Hd = design(f,'iirls','Weights',W);
isstable(Hd)

hfvt = fvtool(Hd, 'Color','w');
legend(hfvt,'IIR Least-Squares','Location', 'NorthWest')
hfvt(2) = fvtool(Hd,'Analysis','phase','Color','white');
legend(hfvt(2),'IIR Least-Squares','Location', 'NorthEast')

%%
% ff = 0:8000;
fs = 16000;
dbPerOct = 3.0; %dB
f_band = [50 4000];


Noct = 2;
offsetdB = 0.377*Noct;


fmid = 10^mean(log10(f_band));
f_edges = [1 f_band(1) fmid f_band(2) fs/2];

f_intpts = (f_band'*2.^(Noct/2*[-1 0 1]))';

a = db2mag( ...
    dbPerOct/log10(2) * log10(f_edges/fmid) );
a([1 end]) = a([2 end-1]);

a_ = db2mag( ...
    dbPerOct/log10(2) * log10(f_intpts/fmid) );
a_([1 end]) = a_(2,:);


a_(f_intpts == f_band) = db2mag(mag2db(a_(f_intpts == f_band)) + [1 -1]'*offsetdB);

f_i_rnd = round(f_intpts);
a_before = a(1)*ones(1,numel(f_edges(1):f_i_rnd(1,1)));
a_int = a_before;
for cf = 1:numel(f_band)
    interpFrqs = f_i_rnd(1,cf)+1:f_i_rnd(end,cf);
    
    a_interp = db2mag( ...
        lagrangepoly( ...
        log10(f_intpts(:,cf)), ...
        mag2db(a_(:,cf)), ...
        log10(interpFrqs))  );
    
    if cf < numel(f_band)
        a_after = db2mag( ...
            dbPerOct/log10(2) * log10( ...
            (f_i_rnd(end,cf)+1:f_i_rnd(1,cf+1)) /fmid) );
    elseif f_i_rnd(end,cf) < f_edges(end)
        a_after = a(end)*ones(1,numel(f_i_rnd(end,end)+1:f_edges(end)));
    else
        a_after=[];
    end
    
    a_int = [a_int a_interp a_after];
end

f = f_edges(1):f_edges(end);

plot(f_edges/1e3,mag2db(a)); hold on;
plot(f/1e3,mag2db(a_int)); hold on;
grid on; grid minor;
set(gca,'xscale','log');
xlim([0.01 10]); hold off;

%%
% [num,den]=iirlpnorm(8,8,f/(fs/2),f/(fs/2),a_int);
[num,den]=iirlpnorm(8,8,f/(fs/2),f/(fs/2),f/(fs/2));

fvtool(num,den);


%%
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
[yy,zz] = meshgrid(linspace(0,3,rtxN)); % Planar Array
% yy = linspace(0,3,rtxN); zz = yy*0+1.5; % Linear Array

rtx = [zeros(numel(yy),1), yy(:), zz(:)];
srx = rtx;

MM=[];PP=[];
tic;
ss=0;
while true %for ss = 1:10
    ss = ss+1;
% r = [1.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
% s = [1.5 1.5 1.5];              % Source position [x y z] (m)
r = rand(1,3)*3;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = rand(1,3)*3;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
% s = [rand(1,2)*3 1.5]; r = [rand(1,2)*3 1.5]; % When using linear array

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

PhaseDifference = mod(unwrap(angle(HF)) - unwrap(angle(HC)) + pi,2*pi)/pi*180-180;
PhaseDifference(end/2+1:end)=[];

MM(:,ss) = MagnitudeC; %.*ff.'; Planar compensation %.*sqrt(ff).' % Linear compensation
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

