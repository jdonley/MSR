
SpeechFilePath = '\+Miscellaneous\+Speech_Files\Female_SA2.WAV';

[x,sig_info.Fs] = audioread(SpeechFilePath);
x = x/max(x(:));

sig_info.Nfft = 12 * 1e-3 * sig_info.Fs; % Number of fft components
sig_info.time_delay = []; % Seconds %If empty the time delay will based on the frame length
sig_info.overlap = 1-1/2;
sig_info.zeropadtime = 0; % miliseconds of zero padding (i.e. maximum time shift before circular convolution)
sig_info.predict_buff = 2; % Length of prediction buffer as a percentage of a full frame prediction (i.e. 10 is %1000 of frame length)
% sig_info.AR_method = 'lpc';
% sig_info.AR_method = 'burg';
sig_info.AR_method = 'cov';
sig_info.window = false;


%%
% falias = 2000;
% [b,a] = cheby1(9,1,falias/(sig_info.Fs/2));
% x = filter(b,a,x);

%%
[ ~, ~, Frames, ~, ~, Frames_orig] = ...
    Broadband_Tools.FFT_custom(x, sig_info);

%%
% for EgTime = (sig_info.Nfft/sig_info.Fs)*3:((1-sig_info.overlap)*sig_info.Nfft/sig_info.Fs):(numel(x)/sig_info.Fs); % seconds
% EgTime = 0.27; % seconds
% EgTime = 0.61; % seconds
% EgTime = 0.90; % seconds
EgTime = 1.24; % seconds


EgFrame = floor(EgTime * sig_info.Fs / ((1-sig_info.overlap) * sig_info.Nfft));
N=1/(1-sig_info.overlap);
% Plot Prediction
fH = figure(1);
linewidth = 1.5;

t = ((1:sig_info.Nfft*3) - sig_info.Nfft*2)/sig_info.Fs*1e3;
plot(t(1:sig_info.Nfft*2),[Frames_orig(EgFrame-N*2,:), Frames_orig(EgFrame-N,:)],'k','linewidth',linewidth);
ax = gca; 
hold on
plot(t(sig_info.Nfft*2+(1:sig_info.Nfft)), Frames_orig(EgFrame,:),':k','linewidth',linewidth);

plot(t(sig_info.Nfft*2+(0:sig_info.Nfft)),[Frames_orig(EgFrame-N,end) Frames((EgFrame-N)*2,:)],'b','linewidth',linewidth);

residual = Frames_orig(EgFrame,:)-Frames((EgFrame-N)*2,:);
plot(t(sig_info.Nfft*2+(1:sig_info.Nfft)),residual,'r','linewidth',linewidth);

% ylim(ax.YLim);
ylim([min(x) max(x)]);
xlim(t([1 end]));

plot(t(sig_info.Nfft*2*[1 1]),ax.YLim,'-k','linewidth',linewidth*2);
xx = (sig_info.Nfft*(1-sig_info.overlap) + [0 sig_info.Nfft*(sig_info.overlap)*[1 1] 0])/sig_info.Fs*1e3;
yy = ax.YLim([1 1 2 2]);
patch(xx,yy,[0 0 0],'edgecolor','none','facealpha',0.25)
hold off

ylabel('Amplitude')
xlabel('Time (ms)')
FontSz = 16;
txtoffset = 1.1;% factor
text(-1*sig_info.Nfft/sig_info.Fs*1e3,max(x)*txtoffset,...
    'Past (Known Signal)',...
    'FontSize',FontSz,'ho','c')

text(0.5*sig_info.Nfft/sig_info.Fs*1e3,max(x)*txtoffset,...
    'Future (Predicted Signal)',...
    'FontSize',FontSz,'ho','c','co','b')

annotation('textarrow',...
    ([(0.20) 0.25]*sig_info.Nfft/sig_info.Fs*1e3  - ax.XLim(1))/diff(ax.XLim)*ax.Position(3) + ax.Position(1),...
    ([min(x)*(1-(txtoffset-1)*3.5)...
    residual(0.25*sig_info.Nfft)] - ax.YLim(1))/diff(ax.YLim)*ax.Position(4) + ax.Position(2),...
    'String','Residual',...
    'FontSize',FontSz,'ho','l','co','r','linewidth',linewidth)


% drawnow;
% pause(1/30);
% end

ax.FontSize = FontSz;
fH.Color = 'w';
fH.Position(3) = fH.Position(4)*2.5;
% drawnow;
% tightfig;