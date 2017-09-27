

%%
% fs = 16000;
% dbPerOct = 3.0; %dB
% f_band = [100 2000];
% 
% 
% Noct = 2;
% offsetdB = 0.377*Noct;
% 
% 
% fmid = 10^mean(log10(f_band));
% f_edges = [1 f_band(1) fmid f_band(2) fs/2];
% 
% f_intpts = (f_band'*2.^(Noct/2*[-1 0 1]))';
% 
% a = db2mag( ...
%     dbPerOct/log10(2) * log10(f_edges/fmid) );
% a([1 end]) = a([2 end-1]);
% 
% a_ = db2mag( ...
%     dbPerOct/log10(2) * log10(f_intpts/fmid) );
% a_([1 end]) = a_(2,:);
% 
% 
% a_(f_intpts == f_band) = db2mag(mag2db(a_(f_intpts == f_band)) + [1 -1]'*offsetdB);
% 
% f_i_rnd = round(f_intpts);
% a_before = a(1)*ones(1,numel(f_edges(1):f_i_rnd(1,1)));
% a_int = a_before;
% for cf = 1:numel(f_band)
%     interpFrqs = f_i_rnd(1,cf)+1:f_i_rnd(end,cf);
%     
%     a_interp = db2mag( ...
%         lagrangepoly( ...
%         log10(f_intpts(:,cf)), ...
%         mag2db(a_(:,cf)), ...
%         log10(interpFrqs))  );
%     
%     if cf < numel(f_band)
%         a_after = db2mag( ...
%             dbPerOct/log10(2) * log10( ...
%             (f_i_rnd(end,cf)+1:f_i_rnd(1,cf+1)) /fmid) );
%     elseif f_i_rnd(end,cf) < f_edges(end)
%         a_after = a(end)*ones(1,numel(f_i_rnd(end,end)+1:f_edges(end)));
%     else
%         a_after=[];
%     end
%     
%     a_int = [a_int a_interp a_after];
% end
% 
% f = f_edges(1):f_edges(end);
% 
% plot(f_edges/1e3,mag2db(a)); hold on;
% plot(f/1e3,mag2db(a_int)); hold on;
% grid on; grid minor;
% set(gca,'xscale','log');
% xlim([0.01 10]); hold off;

%%
fs = 16000;
f_band = [100 2000];
fmid = 10^mean(log10(f_band));
fmid = f_band(end);
res = 100;
F = (0:res:fs/2)/(fs/2);
A = ((0:res:fs/2))/fmid;
H = A .* exp(1j*pi/2);
f = fdesign.arbmagnphase('Nb,Na,F,H',4,1,F,H);
W = [ 0*ones(1,numel(0:res:f_band(1)-1)) ...
      10*ones(1,numel(f_band(1):res:f_band(2) )) ...
      0*ones(1,numel(f_band(2)+1:res:fs/2)) ];
Hd = design(f,'iirls','Weights',W);

isstable(Hd)
Hd.impzlength

imp = Hd.impulse;
imp = imp.Data;

% fvtool(Hd,'polezero')
% fvtool(Hd,'impulse')

% hfvt = fvtool(Hd,'Analysis','freq', 'Fs',16000, 'PhaseUnits','Degrees','Color','w');
% ax = findall(hfvt.Children,'Type','Axes');
% ax.XScale = 'log';

%%
% [num,den]=iirlpnorm(8,8,f/(fs/2),f/(fs/2),a_int);
% [num,den]=iirlpnorm(8,8,f/(fs/2),f/(fs/2),f/(fs/2));
% 
% fvtool(num,den);


%%
% clear; clc;
% close all;


c = 343;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
L = [3 3 3];                % Room dimensions [x y z] (m)
n = 0.2*fs;                 % Number of samples
mtype = 'omnidirectional';  % Type of microphone
mtypeW= 'cardioid';         % Type of microphone
order = 1;                  % -1 equals maximum reflection order
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter
rng shuffle;

%%%
betaA = (1 - [1.0   [1 1 1 1 1]*1.0]).^2;                 % Reverberation time (s)
betaI = (1 - [1.0   [1 1 1 1 1]*1.0]).^2;                 % Reverberation time (s)
%%%
beta1 = (1 - [1.0   [1 1 1 1 1]*1.0]).^2;                 % Reverberation time (s)
% betaW = (1 - [1.0 [1 0 1 1 1]*1.0]).^2;                 % Reverberation time (s)
%%%

rtxN = 24;
[yy,zz] = meshgrid(linspace(0,3,rtxN)); % Planar Array
% yy = linspace(0,3,rtxN); zz = yy*0+1.5; % Linear Array

rtx = [zeros(numel(yy),1), yy(:), zz(:)];
srx = rtx;

MC=[];MF=[];M=[];PP=[];
tic;
ss=0;
while true %for ss = 1:10
    ss = ss+1;
r = [1.0 1.5 1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = [1.5 1.5 1.5];    % Source position [x y z] (m)
% r = rand(1,3)*3;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
% s = rand(1,3)*3;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
% s = [rand(1,2)*3 1.5]; r = [rand(1,2)*3 1.5]; % When using linear array

% hI = rir_generator(c, fs, r, s, L, betaI, n, mtype, order, dim, orientation, hp_filter);
% hA = rir_generator(c, fs, r, s, L, betaA, n, mtype, order, dim, orientation, hp_filter);
% h1 = rir_generator(c, fs, r.*[ 1 1 1], s, L, beta1, n, mtype, order, dim, orientation, hp_filter);
h2 = rir_generator(c, fs, r.*[-1 1 1], s, L, beta1, n, mtype, order, dim, orientation, hp_filter);
% hf = h1(:)+h2(:)-hI(:);
hf = h2(:);

stx = s;              % Source position [x y z] (m)
htx = rir_generator(c, fs, rtx, stx, L, beta1, n, mtype, order, dim, orientation, hp_filter);
rrx = r;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
hrx=[];
for i = 1:size(rtx,1)
    hrx(i,:) = rir_generator(c, fs, rrx, srx(i,:), L, betaA, n, mtype, order, dim, orientation, hp_filter);
end

% htx = imag(hilbert(htx));
hrx = Tools.fconv(hrx.',repmat(imp.',size(hrx,1),1).').';

hc = Tools.fconv(htx.',hrx.');
hc = sum(hc(1:numel(hf),:),2);
% [~,adjV(ss)] = Broadband_Tools.power_norm(hf,hc,fs,[250 1000]);

[b,a] = cheby1(6,0.1,[250 1500]/(fs/2));
% hI_band = filter(b,a,hI);
hf_band = filter(b,a,hf);
hc_band = filter(b,a,hc);
figure(1);
% plot(hI_band); hold on
plot(hf_band); hold on
plot(hc_band); hold on;
% plot(hf_band - hc_band); hold on;
hold off


h = hf-hc;

HF = fft(hf);
HC = fft(hc);
H = fft(h);

ff = linspace(0,fs/2,n/2+1)/1e3;ff(end)=[];

MagnitudeC = abs(HC);
MagnitudeC(end/2+1:end)=[];
MagnitudeF = abs(HF);
MagnitudeF(end/2+1:end)=[];
Magnitude = abs(H);
Magnitude(end/2+1:end)=[];

PhaseDifference = mod(unwrap(angle(HF)) - unwrap(angle(HC)) + pi,2*pi)/pi*180-180;
PhaseDifference(end/2+1:end)=[];

MC(:,ss) = MagnitudeC;
MF(:,ss) = MagnitudeF;
M(:,ss) = Magnitude;
PP(:,ss) = PhaseDifference;
% MM = mean([MM , MagnitudeC.*ff.' ],2);
% PP = mean([PP , PhaseDifference  ],2);

% figure(2);
% subplot(2,1,1);
% plot(ff, mag2db(  MC  ),':k'); hold on;
% plot(ff, mag2db(  mean(MC,2)  ),'-b','linew',1.5); hold on;
% plot(ff, mag2db(  MF  ),':m'); hold on;
% plot(ff, mag2db(  mean(MF,2)  ),'-g','linew',1.5); hold off;
% xlim([0.1 10]); %ylim([-60 0]);
% grid on; grid minor; set(gca,'xscale','log');
% xlabel('Frequency (kHz)');ylabel('Magnitude (dB)');
% 
% subplot(2,1,2);
% plot(ff, PP ,':k'); hold on;
% plot(ff, mean(PP,2) ,'-','co',.8*[1 0 0],'linew',1.5); hold off;
% xlim([0.1 10]); ylim([-180 180]); yticks([-180:45:180]);
% grid on; grid minor; set(gca,'xscale','log');
% xlabel('Frequency (kHz)');ylabel('Phase (\circ)');
hold on
figure(2);
meanMF = mag2db(mean(MF,2));
% plot(ff, mag2db(  MF        ) - meanMF  ,':k'); hold on;
% plot(ff, mag2db(  M         ) - meanMF  ,':m'); hold on;
plot(ff, mag2db(  mean(MF,2)) - meanMF  ,'-b','linew',1.5); hold on;
plot(ff, mag2db(  mean(M,2) ) - meanMF  ,'-r','linew',1.5); hold off;
xlim([0.1 10]); %ylim([-60 0]);
grid on; grid minor; set(gca,'xscale','log');
xlabel('Frequency (kHz)');ylabel('Magnitude (dB)');

drawnow;
fprintf('%d\n',ss);
end
toc
