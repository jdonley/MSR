
SpeechFilePath = '\+Miscellaneous\+Speech_Files\Female_SA2.WAV';

[x,sig_info.Fs] = audioread(SpeechFilePath);

sig_info.Nfft = 12 * 1e-3 * sig_info.Fs; % Number of fft components
sig_info.time_delay = []; % Seconds %If empty the time delay will based on the frame length
sig_info.overlap = 1-1/2^5;
sig_info.zeropadtime = 0; % miliseconds of zero padding (i.e. maximum time shift before circular convolution)
sig_info.predict_buff = 2; % Length of prediction buffer as a percentage of a full frame prediction (i.e. 10 is %1000 of frame length)
sig_info.AR_method = 'burg';
sig_info.window = false;


%%
% falias = 2000;
% [b,a] = cheby1(9,1,falias/(sig_info.Fs/2));
% x = filter(b,a,x);

%%
[ ~, ~, Frames, ~, ~, Frames_orig] = ...
    Broadband_Tools.FFT_custom(x, sig_info);

%%
for EgTime = (sig_info.Nfft/sig_info.Fs)*3:((1-sig_info.overlap)*sig_info.Nfft/sig_info.Fs):(numel(x)/sig_info.Fs); % seconds
% EgTime = 0.27; % seconds
% EgTime = 0.61; % seconds
% EgTime = 0.90; % seconds


EgFrame = floor(EgTime * sig_info.Fs / ((1-sig_info.overlap) * sig_info.Nfft));
N=1/(1-sig_info.overlap);
%% Plot Prediction
linewidth = 1.5;

t = ((1:sig_info.Nfft*3) - sig_info.Nfft*2)/sig_info.Fs*1e3;
plot(t(1:sig_info.Nfft*2),[Frames_orig(EgFrame-N*2,:), Frames_orig(EgFrame-N,:)],'k','linewidth',linewidth);
ax = gca; 
hold on
plot(t(sig_info.Nfft*2+(1:sig_info.Nfft)), Frames_orig(EgFrame,:),':k','linewidth',linewidth);

plot(t(sig_info.Nfft*2+(0:sig_info.Nfft)),[Frames_orig(EgFrame-N,end) Frames((EgFrame-N)*2,:)],'r','linewidth',linewidth);

% ylim(ax.YLim);
ylim([min(x) max(x)]);
xlim(t([1 end]));

plot(t(sig_info.Nfft*2*[1 1]),ax.YLim,'-k','linewidth',linewidth);
xx = (sig_info.Nfft*(1-sig_info.overlap) + [0 sig_info.Nfft*(sig_info.overlap)*[1 1] 0])/sig_info.Fs*1e3;
yy = ax.YLim([1 1 2 2]);
patch(xx,yy,[0 0 0],'edgecolor','none','facealpha',0.1)
hold off

ylabel('Amplitude')
xlabel('Time (ms)')

drawnow;
pause(1/30);
end