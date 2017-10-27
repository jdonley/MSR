
 clear; 

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
c = 343;
rtxN = 60;
startL = 0;
endL = 3;
linePos = (startL + endL/rtxN/2) : endL/rtxN : endL*(1 - 1/rtxN/2);
[yy,zz] = meshgrid( linePos ); % Planar Array
d = mean(mean([diff(zz,[],1) diff(yy,[],2).'] ));

f_lo = c / (2*endL);
f_hi = c / (2*d);


nb = 4;
na = 1;

s_up_fact = 3;

fs = 16000*s_up_fact;
f_band = [25 f_hi];%*s_up_fact;
% f_band = round([f_lo f_hi]);
f_filtlow = 10*s_up_fact;
% fmid = 10^mean(log10(f_band));

% [bc,ac]=cheby1(6,1,f_band(2)/(fs/2));
% impc = impz(bc,ac);
% IMPC = fft(impc,1024*2);
% IMPC(end/2+1:end) = [];
% Ac = -unwrap(angle(IMPC));

% res = 20;
% F = [0 ...
%     (res:res:fs/2)/(fs/2)];
F = [0 ...
    logspace(log10(f_filtlow),log10(fs/2),1023)/(fs/2)];  F(end)=1;
w = [1 ... DC component
     1*ones(1,numel(F(F< f_band(1)/(fs/2)) ) - 1) ...
     1*ones(1,numel(F(F>=f_band(1)/(fs/2) & F<=f_band(2)/(fs/2)))) ...
     0*ones(1,numel(F(F> f_band(2)/(fs/2)))) ...
     ];
A = F*fs/2;%[0 ... DC component
    %res:res:fs/2];
P = [0 ... DC component
	...Ac(2:end).' + ...
    ones(1,numel(A)-1)*pi/2];
% ff = linspace(F(1),F(end),nfft);
% Pbegin = P((ff<f_band(1)/(fs/2)) );
% Pend = P((ff>f_band(2)/(fs/2)));
% wp = [ ... DC component
%      Pbegin(end)+Pbegin*0 ...
%      P((ff>=f_band(1)/(fs/2) & ff<=f_band(2)/(fs/2))) ...
%      Pend(1)+Pend*0 ...
%      ];
% P = wp; P(1) = 0;

H = A .* exp(1j*P);
f = fdesign.arbmagnphase('Nb,Na,F,H',nb,na,F,H);


Hd = design(f,'iirls','Weights',w);

imp_ = Hd.impulse;
imp_ = imp_.Data;

if isstable(Hd), ImpSt='true';else,ImpSt='false';end
fprintf('WFS/SDM IIR(LS) pre-filter is stable: %s\n',ImpSt);
fprintf('WFS/SDM IIR(LS) pre-filter length: %d\n',numel(imp_));



% fvtool(Hd,'polezero')
% fvtool(Hd,'impulse')

% hfvt = fvtool(Hd,'Analysis','freq', 'Fs',16000, 'PhaseUnits','Degrees','Color','w');
% ax = findall(hfvt.Children,'Type','Axes');
% ax.XScale = 'log';

%%% 
nfft = max(nextpow2(numel(H)),4096);
ff = linspace(F(1),F(end),nfft);
aa = interp1(F,A,ff);
pp = interp1(F,P,ff);
HH = aa.*exp(1i*pp);
WW = interp1(F,w,ff);
% WW(WW~=0) = tukeywin(nnz(WW),0.1);

WW = sqrt(WW);
OM = exp(-1i*(0:nb)' * ff*pi);
Dva =  (OM(2:na+1,:).') .* HH.';
Dvb = -(OM(1:nb+1,:).');
% D=[Dva Dvb].*(WW.'*ones(1,na+nb+1));
% 
% R  = real(D'*D);
% Vd = real(D'*(-HH.*WW).');
% 
% th = R\Vd;
% 
% D_ = diag(WW)*[Dvb Dva];
% HH_ = diag(WW)*(-HH).';
% th = real(D_'*D_) \ real(D_' * HH_);

OM = exp(-1i*(0:nb)' * ff*pi);
Dvb = -(OM(1:nb+1,:).');
Dva =  (OM(2:na+1,:).') .* HH.';

D=[Dvb Dva];
W = diag(WW);
th = real( D' * W * D ) \ real( D' * W * (-HH).' );

b = th(1:nb+1).';
a = [1 th(nb+2:end).'];

a = polystab(a);

imp = impz(b,a);
impIDEAL = (ifft([HH conj(HH(end:-1:2))])).'; % But contains delay
% imp(mag2db(abs(imp/max(imp)))<-120) = [];

if isstable(b,a), ImpSt='true';else,ImpSt='false';end
fprintf('WFS/SDM IIR(LS) pre-filter is stable: %s\n',ImpSt);
fprintf('WFS/SDM IIR(LS) pre-filter length: %d\n',numel(imp));


PRE_b = b;
PRE_a = a;
0;

% imp = Tools.fconv(imp,impc);

% y = zeros(1,16000);
% y(end/2)=1;

% yy = Tools.fconv(y.',imp);
% yy2 = filter(b,a,y).';
% yy3 = Tools.fconv(y.',imp);
% 
% Y   = fft(y);
% YY2 = fft(yy2(1:numel(y)).');
% YY3 = fft(yy3(1:numel(y)).');
% Y(end/2:end) = [];
% YY2(end/2:end) = [];
% YY3(end/2:end) = [];

IMP = fft(imp,4096);
IMP(end/2+1:end) = [];
% IMP_ = fft(imp_,1024);
% IMP_(end/2+1:end) = [];
frqs = linspace(0,fs/2,numel(IMP))/1e3;

close all;
fH = figure(1);
fH.Color = 'w';
yyaxis left;
ax = gca;
blueColor = [0.0 0.3 0.7];
plot(frqs,mag2db(abs(IMP)),'.','color',blueColor,'linew',1.5); hold on;
% plot(frqs,mag2db(abs(IMP_)),'.','color',magColor/2,'linew',1.5); hold on;
plot(ff*fs/2/1e3,WW.^2 * 99+0.5,'color','k','linew',1.5); hold on;
hold off;

text(1.9,20,'$k_{\mathrm{u}}$','Interpreter','latex','fontsize',12);
ar = annotation(gcf,'arrow', [0.58 0.615],[0.36 0.26]);
ar.HeadStyle = 'deltoid';


ax.YAxis(1).Label.String = {'Magnitude (dB)';  'or  LS Weight (\%)'};
ax.YAxis(1).Label.Interpreter = 'latex';
ax.YAxis(1).Color = blueColor;
ax.YAxis(1).MinorTick = 'on';
ax.YAxis(1).TickDirection = 'both';
ylim([0 100]); 

yyaxis right;
ax = gca;
redColor = [0.8 0.1 0.1];
plot(frqs,(mod(unwrap(angle(IMP))+pi,2*pi)-pi)/pi*180,'.','color',redColor,'linew',1.5); hold on
% plot(frqs,(mod(unwrap(angle(IMP_))+pi,2*pi)-pi)/pi*180,'.','color',phaseColor/2,'linew',1.5); hold on
% plot(frqs,unwrap(mod(unwrap(angle(Y)) - unwrap(angle(YY2)+pi),2*pi)-pi)/pi*180); hold on
% plot(frqs,unwrap(mod(unwrap(angle(Y)) - unwrap(angle(YY3)+pi),2*pi)-pi)/pi*180); hold on
% plot(f_band/1e3,[90 90],'-k','linew',1.5);  hold on
plot(f_band/1e3,[91 91],':k','linew',0.5);  hold on
plot(f_band/1e3,[89 89],':k','linew',0.5);  hold on
hold off;
ax.YAxis(end).Label.String = 'Phase (${}^\circ$)';
ax.YAxis(end).Label.Interpreter = 'latex';
ax.YAxis(end).Color = redColor;
ax.YAxis(end).MinorTick = 'on';
ax.YAxis(end).TickDirection = 'both';
ax.XScale = 'log';
ax.XAxis.TickDirection = 'both';
ax.XAxis(1).Label.String = 'Frequency (kHz)';
ax.XAxis(1).Label.Interpreter = 'latex';
ylim([60 120]); 
xlim([0.1 fs/2/1e3])
ax.XTick = [0.1, 1, 3.43, 10, 24];
ax.XTickLabel = {'0.1','1','3.43','10','24'};
grid off; grid on; grid minor;
fH.Units = 'centimeters';
fH.Position(3:4) = [12 5];
tightfig;

%%%
fH2 = figure(2);
fH2.Color = 'w';

yyaxis right;
ax = gca;
stem((0:numel(imp)-1)/fs*1e3,...
    mag2db(abs(imp)),...
    ':x','color', redColor);
ylim( max(mag2db(abs(imp)))*[-1 1]*1.1 );
ax = gca;
ax.YAxis(end).Label.String = 'Magnitude (dB)';
ax.YAxis(end).Label.Interpreter = 'latex';
ax.YAxis(end).Color = redColor;
ax.YAxis(end).MinorTick = 'on';
ax.YAxis(end).TickDirection = 'both';

yyaxis left;
ax = gca;
stem((0:numel(imp)-1)/fs*1e3,...
    ((imp)),...
    'color', blueColor);
ylim( max(abs(imp))*[-1 1]*1.1 );
ax = gca;
ax.YAxis(1).Label.String = 'Amplitude';
ax.YAxis(1).Label.Interpreter = 'latex';
ax.YAxis(1).Color = blueColor;
ax.YAxis(1).MinorTick = 'on';
ax.YAxis(1).TickDirection = 'both';

ax.XLim = ax.XLim + [-1 1]*diff(ax.XLim)*0.03;
ax.XAxis.TickDirection = 'both';
ax.XAxis.MinorTick = 'on';
ax.XAxis.Label.String = 'Time (ms)';
ax.XAxis.Label.Interpreter = 'latex';
grid off; grid on; 
fH2.Units = 'centimeters';
fH2.Position(3:4) = [12 5];
tightfig;
%%
% [num,den]=iirlpnorm(8,8,f/(fs/2),f/(fs/2),a_int);
% [num,den]=iirlpnorm(8,8,f/(fs/2),f/(fs/2),f/(fs/2));
%
% fvtool(num,den);


%%
 % clc;
% close all;


% c = 343;                    % Sound velocity (m/s)
% fs = 16000;                 % Sample frequency (samples/s)
L = [3 3 3];                % Room dimensions [x y z] (m)
n = 0.1*fs;                 % Number of samples
mtype = 'omnidirectional';  % Type of microphone
mtypeW= 'cardioid';         % Type of microphone
order = 2;                  % -1 equals maximum reflection order
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter
rng(sum('ICASSP2018'));

%%%
betaW     = (1 - [1.0   [1 1 1 1 1]*1.0]).^2;   %anechoic              % Reverberation time (s)
%%%
beta(1,:) = (1 - [1.0   [1 1 1 1 1]*1.0]).^2;                 % Reverberation time (s)

beta(2,:) = (1 - [1.0   [0 1 1 1 1]*1.0]).^2;                 % Reverberation time (s)
beta(3,:) = (1 - [1.0   [0 0 1 1 1]*1.0]).^2;                 % Reverberation time (s)
beta(4,:) = (1 - [1.0   [0 0 0 1 1]*1.0]).^2;                 % Reverberation time (s)
beta(5,:) = (1 - [1.0   [0 0 0 0 1]*1.0]).^2;                 % Reverberation time (s)
beta(6,:) = (1 - [1.0   [0 0 0 0 0]*1.0]).^2;                 % Reverberation time (s)
%%%

% rtxN = 60;
% % linePos = linspace(0,3,rtxN);
% startL = 0;
% endL = 3;
% linePos = (startL + endL/rtxN/2) : endL/rtxN : endL*(1 - 1/rtxN/2);
% [yy,zz] = meshgrid( linePos ); % Planar Array
% % yy = linspace(0,3,rtxN); zz = yy*0+1.5; % Linear Array
% 
% d = mean(mean([diff(zz,[],1) diff(yy,[],2).'] ));

rtx = [zeros(numel(yy),1), yy(:), zz(:)]; % monopole
% rtx = reshape( ...
%        [zeros(numel(yy),1)+0.025.*[-1 1], ...
%        [yy(:) yy(:)], ...
%        [zz(:) zz(:)]], ...
%        [], 3); % dipole

srx = rtx;

imgs = 6;
% res = 20;
% [XX,YY] = meshgrid(linspace(0,3,3*res));


%%% taper window to limit diffraction
winPerc = 20;
[Wx,Wy] = meshgrid( tukeywin(rtxN,winPerc/100) );
DiffracWin = Wx .* Wy;
%%%


% [b,a] = cheby1(6,1,[100 2000]/(fs/2));




MC=[];MF=[];M=[];MI=[];PP=[];
% hh1=zeros((3*res)^2,n);
% hh2=hh1;
tic;
ss=0;
while ss < 10%200 %ss<1 %for ss = 1:10
    ss = ss+1;
% for ss = 1%:(3*res)^2
%     [x_,y_] = ind2sub(size(XX),ss);
%     
%     x = XX(x_,y_); y = YY(x_,y_);
% x=1.0; 
% y=1.5;
    r = rand(1,3).*[2.5 3 3] + [0.5 0 0];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
    s = rand(1,3).*[2.5 3 3] + [0.5 0 0];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
    % s = [rand(1,2)*3 1.5]; r = [rand(1,2)*3 1.5]; % When using linear array
%     r  = [ 0.5   1.5   1.5];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
%     s  = [ 1.5   1.5   1.5];    % Source position [x y z] (m)
    

    for img = imgs%1:6
        
        %%% Mic transfer functions
        stx = s;              % Source position [x y z] (m)
        htx = rir_generator(c, fs, rtx, stx, L, beta(img,:), n, mtype, order-1, dim, orientation, hp_filter);
        htxLR = htx - ... % Last reflection
            rir_generator(c, fs, rtx, stx, L, beta(img,:), n, mtype, order-2, dim, orientation, hp_filter) ;
        %%%

        %%% Ground truth reflections
        h1 = rir_generator(c, fs, r, s, L, [1 beta(img,2:end)], n, mtype, order, dim, orientation, hp_filter);
        h2 = rir_generator(c, fs, r, s, L, [0 beta(img,2:end)], n, mtype, order, dim, orientation, hp_filter);
        h1_ = h1-h2;
        hf = h1_(:);
        %%%


        %%% Loudspeaker transfer functions
        rrx = r;    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
        hrx=[]; hrxLR=[];
        for i = 1:size(rtx,1)
            hrx(i,:) = rir_generator(c, fs, rrx, srx(i,:), L, beta(img,:), n, mtype, order-1, dim, orientation, hp_filter);
            hrxLR(i,:) = rir_generator(c, fs, rrx, srx(i,:), L, beta(img,:), n, mtype, order-2, dim, orientation, hp_filter);
        end
        
        %%% Determine cardioid pattern multiplier
        th = atan( ...
            sqrt( (srx(:,2) - rrx(2)).^2 + (srx(:,3) - rrx(3)).^2 ) ...
            ./ (srx(:,1) - rrx(1)) );
%         alph = 0.5; % cardioid
%         A = alph + (1-alph)*cos(th);  
%         alph = realmax; % monopole
%         alph = 2/1.0; % sub-cardioid
%         alph = 1/1.0; % cardioid
%         alph = 1/1.7; % super-cardioid
%         alph = 1/3.0; % hyper-cardioid     
        alph = 0/1.0; % dipole   
        A =  (alph + cos(th)) ...
            /(abs(alph) + 1);
        %%% Apply directional pattern to microphones and loudspeakers
        htxDi = htx .* A;
        htxLRDi = htxLR .* A;          
%         hrx = hrx .* A;
%         hrxLR = hrxLR .* A;  

%         htx = htx(1:end/2,:) + htx(end/2+1:end,:);
%         htxLR = htxLR .* A;          
%         hrx = hrx .* A;
%         hrxLR = hrxLR .* A;    
        %%%
        
        %%% Apply WFS/SDM pre-filter
        hrx = Tools.fconv(hrx.',repmat(imp.',size(hrx,1),1).').';
        hrxLR = Tools.fconv(hrxLR.',repmat(imp.',size(hrxLR,1),1).').';
        
%         hrxIdeal = Tools.fconv(hrx.',repmat(impIDEAL.',size(hrx,1),1).').';
%         hrxLRIdeal = Tools.fconv(hrxLR.',repmat(impIDEAL.',size(hrxLR,1),1).').';
        %%%
        
        %%% Cancellation signal minus last reflection
        %%% ONE (previously suggested pre-filter)
%         hcId = Tools.fconv(htx.',hrxIdeal.');        
%         hcLRdirectId = Tools.fconv(htxLR.',hrxLRIdeal.');
%         hcLRId = Tools.fconv(htxLR.',hrxIdeal.');
%         hcLId = (hcLRId - hcLRdirectId);
%         %         hc = hc .* repmat(DiffracWin(:).',size(hc,1),1);
% %         hcL = hcL .* repmat(DiffracWin(:).',size(hcL,1),1);
%         hcId = sum(hcId(1:numel(hf),:),2) / rtxN^2 / pi;
%         hcLId = sum(hcLId(1:numel(hf),:),2) / rtxN^2 / pi;
        
        %%% TWO
        hc = Tools.fconv(htx.',hrx.');        
        hcLRdirect = Tools.fconv(htxLR.',hrxLR.');
        hcLR = Tools.fconv(htxLR.',hrx.');
        hcL = (hcLR - hcLRdirect);
        %         hc = hc .* repmat(DiffracWin(:).',size(hc,1),1);
%         hcL = hcL .* repmat(DiffracWin(:).',size(hcL,1),1);
        hc = sum(hc(1:numel(hf),:),2) / rtxN^2 / pi;
        hcL = sum(hcL(1:numel(hf),:),2) / rtxN^2 / pi;
        
        %%% THREE
        hcDi = Tools.fconv(htxDi.',hrx.');        
        hcLRdirectDi = Tools.fconv(htxLRDi.',hrxLR.');
        hcLRDi = Tools.fconv(htxLRDi.',hrx.');
        hcLDi = (hcLRDi - hcLRdirectDi);
    %         hc = hc .* repmat(DiffracWin(:).',size(hc,1),1);
%         hcL = hcL .* repmat(DiffracWin(:).',size(hcL,1),1);
        hcDi = sum(hcDi(1:numel(hf),:),2) / rtxN^2 / pi;
        hcLDi = sum(hcLDi(1:numel(hf),:),2) / rtxN^2 / pi;    
        
        

        
        
%         hcId = hcId - hcLId;
        hc = hc - hcL;
        hcDi = hcDi - hcLDi;
        %%% 
        
%         % hI_band = filter(b,a,hI);
%         hf_band = filter(b,a,hf);
% %         h1_band = filter(b,a,(h1 ));
% %         h2_band = filter(b,a,hctest_reflect);
%         hc_band = filter(b,a,hc);
% %         hcL_band = filter(b,a,hcL);
%         
%         
%         figure(2); 
%         % plot(hI_band); hold on
%         plot(hf_band); hold on
% %         plot(hf,'k'); hold on
% %         plot(h2,'r'); hold on
%         plot(hc_band); hold on;
% %         plot(hcL_band); hold on;
% %         plot(h1_band); hold on;
% %         plot(h2_band); hold on;
% %         plot(hf_band - hc_band); hold on;
%         hold off;grid on;
% %         pow2db(sum((hf_band).^2)) - pow2db(sum((hf_band-hc_band).^2))
%         0;
        
% hh1(ss,:) = hf_band;
% hh2(ss,:) = hc_band;
        
%         hId = hf-hcId;
        h = hf-hc;
        hDi = hf-hcDi;
        
        HF = fft(hf);
%         HC = fft(hc);
%         HId = fft(hId);
        H = fft(h);
        HDi = fft(hDi);
%         H2 = fft(h2);
%         HI = fft(hI);
        
        ff = linspace(0,fs/2,n/2+1)/1e3;ff(end)=[];
        
%         MagnitudeC = abs(HC);
%         MagnitudeC(end/2+1:end)=[];
        MagnitudeF = abs(HF);
        MagnitudeF(end/2+1:end)=[];
%         MagnitudeId = abs(HId);
%         MagnitudeId(end/2+1:end)=[];
        Magnitude = abs(H);
        Magnitude(end/2+1:end)=[];
        MagnitudeDi = abs(HDi);
        MagnitudeDi(end/2+1:end)=[];
%         Magnitude2 = abs(H2);
%         Magnitude2(end/2+1:end)=[];
        
%         PhaseDifference = mod(unwrap(angle(HF)) - unwrap(angle(HC)) + pi,2*pi)/pi*180-180;
%         PhaseDifference(end/2+1:end)=[];
        
%         MC(:,ss,img) = MagnitudeC;
        MF(:,ss,img) = MagnitudeF;
%         MId(:,ss,img) = MagnitudeId;
        M(:,ss,img) = Magnitude;
        MDi(:,ss,img) = MagnitudeDi;
%         M2(:,ss,img) = Magnitude2;
%         PP(:,ss,img) = PhaseDifference;
        % MM = mean([MM , MagnitudeC.*ff.' ],2);
        % PP = mean([PP , PhaseDifference  ],2);
        
    end
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
    
    fHres = figure(222);
    fHres.Color = 'w';
    meanMF = mean(mean( mag2db(MF),2));
    % plot(ff, mag2db(  MF        ) - meanMF  ,':k'); hold on;
    % plot(ff, mag2db(  M         ) - meanMF  ,':m'); hold on;
    if ss > 1
        CIs = Tools.confidence_intervals(db2mag(mag2db(  MF(:,:,img) ).' - meanMF(:,:,img).'),95,true);
        CIs = mag2db(exp(CIs));
        CI = CIs + mean( mag2db(  MF(:,:,img)),2 ) - meanMF(:,:,img);
        plot(ff, CI  ,'-','color',[0 0 0 0.2],'linew',1.5); hold on;
    end
    pl0 = plot(ff, mean(mag2db(  MF(:,:,img(1))),2 ) - meanMF(:,:,img(1))  ,'-','color',[0 0 0],'linew',1.5); hold on;
    set(gca,'ColorOrderIndex',1);
    for img = imgs%1:6 
        
%         if ss > 1
%             CIs = Tools.confidence_intervals(db2mag(mag2db(  MId(:,:,img) ).' - meanMF(:,:,img).'),95,true);
%             CIs = mag2db(exp(CIs));
%             CI = CIs + mag2db(  mean( MId(:,:,img),2) ) - meanMF(:,:,img);
%             plot(ff, CI  ,'-','color',[1.0 0 1.0 0.2],'linew',1.5); hold on; 
%         end
%         plot(ff, mag2db(  mean( MId(:,:,img),2) ) - meanMF(:,:,img)  ,'-m','linew',1.5); hold on;
        
        if ss > 1
            CIs = Tools.confidence_intervals(db2mag(mag2db(  M(:,:,img) ).' - meanMF(:,:,img).'),95,true);
            CIs = mag2db(exp(CIs));
            CI = CIs + mean( mag2db(  M(:,:,img)),2 ) - meanMF(:,:,img);
            plot(ff, CI  ,'-','color',[redColor 0.2],'linew',1.5); hold on; 
        end
        pl1 = plot(ff, mean( mag2db(  M(:,:,img)),2 ) - meanMF(:,:,img)  ,'-','color',redColor,'linew',1.5); hold on;
        
        if ss > 1
            CIs = Tools.confidence_intervals(db2mag(mag2db(  MDi(:,:,img) ).' - meanMF(:,:,img).'),95,true);
            CIs = mag2db(exp(CIs));
            CI = CIs + mean(  mag2db( MDi(:,:,img)),2 ) - meanMF(:,:,img);
            plot(ff, CI  ,'-','color',[blueColor 0.2],'linew',1.5); hold on; 
        end
        pl2 = plot(ff, mean(  mag2db( MDi(:,:,img)),2 ) - meanMF(:,:,img)  ,'-','color',blueColor,'linew',1.5); hold on;
        
    end
    
    hold off;
    xlim([0.1 fs/2/1e3]); 
    ylim([-20 24]);
    
    lH = legend([pl0, pl1, pl2], ...
        {'Inactive'; ...
        %         'Active Wall On & 1 Reflection'; ...
        %         'Active Wall On & 2 Reflections'; ...
        %         'Active Wall On & 3 Reflections'; ...
        %         'Active Wall On & 4 Reflections'; ...
        %         'Active Wall On & 5 Reflections'; ...
        'Active - WFR: Proposed WLS'; ...
        'Active - Proposed FOD'}, ...
        'Location','northwest', ...
        'interpreter','latex');
    ax = gca;
    ax.XGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    ax.YAxis.Label.String = 'Magnitude (dB)';
    ax.YAxis.Label.Interpreter = 'latex';
    ax.YAxis.MinorTick = 'off';
    ax.YAxis.TickDirection = 'both';
    ax.YTick = [-20:5:24];
    ax.YTickLabel = num2cell([-20:5:24]);
    ax.XScale = 'log';
    ax.XAxis.TickDirection = 'both';
    ax.XAxis.Label.String = 'Frequency (kHz)';
    ax.XAxis.Label.Interpreter = 'latex';
    ax.XTick = [0.1, 1, 3.43, 10, 24];
    ax.XTickLabel = {'0.1','1','3.43','10','24'};
    
    fHres.Units = 'centimeters';
    fHres.Position(3:4) = [12 5];

% figure(111); scatter(x,y,'ok'); hold on;
% xlim([0 3]); ylim([0 3]);



drawnow;

disp(ss);

end
toc;

ax = gca;
ax.ColorOrder = [0 0 0];
yyaxis right;
ax = gca;
ax.YLabel.String = 'Absorption Coefficient';
ax.YLabel.Interpreter = 'latex';
ax.YTick = ([-20:5:20] + 20)/(24+20);
ax.YTickLabel = mat2cell(round( 1-sqrt(db2mag([-20:5:20])) ,2 ) ,1,ones(1,9));
ax.TickDir = 'both';
ax.YMinorTick = 'off';

fHres.Position(3:4) = [12 5];
tightfig;

text(5.5,-11.5,'$k_{\mathrm{u}}$','Interpreter','latex','fontsize',12,'ho','c');
ar = annotation(gcf,'arrow', [0.725 0.685],[0.3 0.2]);
ar.HeadStyle = 'deltoid';

    lH = legend([pl0, pl1, pl2], ...
        {'Inactive'; ...
        %         'Active Wall On & 1 Reflection'; ...
        %         'Active Wall On & 2 Reflections'; ...
        %         'Active Wall On & 3 Reflections'; ...
        %         'Active Wall On & 4 Reflections'; ...
        %         'Active Wall On & 5 Reflections'; ...
        'Active - WFR: Proposed WLS'; ...
        'Active - Proposed FOD'}, ...
        'Location','northwest', ...
        'interpreter','latex');
%%
% hh(:,:,:,1) = reshape(hh1,3*res,3*res,[]);
% hh(:,:,:,2) = reshape(hh2,3*res,3*res,[]);
% 
% [~,IC] = max(abs(hh(ceil(size(XX,1)/2),ceil(size(XX,2)/2),:,1))); % spatial-centre time-index
% FIELDERROR = diff(hh(:,:,IC,:),[],4);
% imagesc(FIELDERROR);

%%
