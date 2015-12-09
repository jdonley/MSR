close all;
f1 = 20; %Hz
f2 = 20000; %Hz
fs = 48000; %Hz
c = 343; %metres/sec
len = 30; %sec
tail = 5; %sec

y = Tools.synthSweep(len+tail,fs,f1,f2,tail*fs);
invFFT = Tools.invSweepFFT(y,f1,f2,fs);
y = [zeros(tail*fs,1); y(:)];
%audiowrite('sinesweep_20Hz-20kHz_30sec_5secTail.WAV',y,fs);

%[y2, fs2] = audioread('result.WAV');
hold on;
delay_avg = []; eq=[];

for dist = 1:2
    for s = 1:2
        for e_ = 1
%  for spkr = 2
%  for rec = 1
    
%     data = audioread(['+Results\+Test_Recordings\+23_11_2015\sinesweep_20Hz-20kHz_30sec_5secTail_spkr' num2str(spkr) '_dist' num2str(dist) '_rec' num2str(rec) '.wav']);
    switch s
        case 1
            Shift = 'no';
        case 2
            Shift = 'with';
    end
    EQ = 'no';
    data = audioread(['+Results\+Test_Recordings\+24_11_2015\distance' num2str(dist) '_' EQ 'EQ_' Shift 'Shift.wav']);
    
    y2 = data(1:(length(y)-(tail-1)*fs));
    
%     f1_=f1;
%     X = fft(y2);
%     % Frequencies
% M = length(X);
% NumPts = M/2 + 1;
% freqs = linspace(0, fs/2, NumPts);
% 
%     W = FrequencyResponse((freqs>=f1_) & (freqs<=f2));
%      I=find(freqs>1000,1,'first');
%      W =  (W ./ W(I));
%      %W = smooth(W,1000);
%             
%      W = -mag2db(W);
%     
% %     plot(freqs((freqs>=f1_) & (freqs<=f2)),W);
% %     grid on;    
% %     set(gca,'YScale','lin');
% %     set(gca,'XScale','log');
%     
% % Weight Levels
% W = [linspace(0, W(1), length(freqs(freqs<f1_))), ...
%      W(:)', ...
%      linspace(W(end), 0, length(freqs(freqs>f2)))];
%  
%  
%  
% 
% % Apply magnitude weighting
% X(1:NumPts) = X(1:NumPts) .* db2mag(W(:));
% 
% % Apply conjugation for negative frequency side of spectrum
% X(NumPts+1:M) = conj(X(M/2:-1:2));
% 
% % Time Domain Transform
% data_ = ifft(X); % Inverse Fast Fourier Transform
% 
% % prepare output vector y
% data_ = real(data_(1:M));
% 
% % remove DC
% data_ = data_(:) - mean(data_);
%     y2=data_;
    
    %y2 = y2 ./ max(abs(y2(:)));
    
    [ir,irinv] = Tools.extractIR(y2,invFFT);
    
    ir_ = ir((tail*fs):end);
    
    irfull = [irinv; ir];
    Zr = fft(irfull);
    fr = abs(Zr( 1:((end/2)-1) ));
    
    %eq((spkr-1)*5+rec,:) = fr;
    
    f = linspace(0,fs/2,length(Zr)/2);f(end)=[];
    
    figure(1); hold on;
    plot(f,(mag2db(fr)));
    set(gca,'XScale','log');
    xlim([20,20000]);
    ylim([-80,0]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    grid on;
    
    figure(2); hold on;
    plot((1:length(ir_))/fs*c, ir_); 
    [pkv,pk]=max(ir_(:));    
    plot(pk/fs*c,pkv,'or');
    xlim([0, 3.0]);
    xlabel('Distance (metres)');
    grid on; grid minor;
    delay_avg(dist,s) = pk;
%  end
%  end
end
end
end
hold off;
delay_ = diff(mean(delay_avg,2));
distance_ = delay_/fs * c; % metres
figure(2);
title(['Distance moved: ' num2str(round(distance_*100,2)) ' centimetres']);

%%
% L = 24;
% b = setup.Loudspeaker_Locations(L,2);
% c = 0.6;
% A = 270-setup.Loudspeaker_Locations(L,1)/pi*180;
% a = sqrt(b^2 + c^2 - 2*b*c*cosd(A));
% 
% B = asind( sind(A) * b / a );
% 
% maxAngle = 90-B


% c = 343;
% phiL = 90/180*pi;
% R = 0.9;
% 
% L = 1:24;
% 
% f = (L-1)*c / (2*R*phiL);
% 
% plot(f,L)
% 
% f = f(1):f(end);
% 
% L = ceil( ((phiL/pi*180) / 360) * 2 * ( ceil( (f/c *2*pi) * R )) + 1 );
% 
% hold on;
% 
% plot(f,L)
% 
% hold off






















%%
% load('LUT_Weight_vs_Frequency_15deg_512f_32w.mat')
% load('LUT_Weight_vs_Frequency_0deg_512f_32w_zones_swapped.mat')
% load('LUT_Weight_vs_Frequency_15deg_512f_32w_compact.mat')
% load('LUT_Weight_vs_Frequency_0deg_512f_32w_zones_swapped_compact.mat')
% load('Database_90degB_-90degQ_0deg_512f_32w.mat')
% load('Database_-90degB_90degQ_-60deg_512f_32w.mat')

%LUT_MagDiff = (Bright_Sample__Weight_Vs_Frequency ./ Quiet_Sample__Weight_Vs_Frequency);
% LUT_MagDiff = mag2db(Acoustic_Contrast__Weight_Vs_Frequency);
% 
% Number_of_Weights = 256;
% len = 4096;
% Fs=16000;
% noise_freqs = linspace(0, Fs/2, len/2 + 1);
% noise_freqs = noise_freqs(noise_freqs>=min(Frequencies) & noise_freqs<=max(Frequencies));
% weights = [0 logspace(log10( 1e-2), log10(  1e4), Number_of_Weights - 1) ];
% 
% LUT_MagDiff_interp = interp2(Frequencies,Weights,LUT_MagDiff,noise_freqs',weights,'spline');
% 
% [val,I]=max(LUT_MagDiff_interp);
% noise_weights = weights(I);
% 
% 
% figure(1);
% surf(noise_freqs,weights, LUT_MagDiff_interp,'linestyle','none');
% view(2);
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% set(gca,'TickDir','out');
% ylabel('Quiet Zone Weight');
% xlabel('Frequency (Hz)');
% c=colorbar;
% ylabel(c,'Zone Contrast (dB)')
% %caxis([0 1]);
% hold on;
% 
%  plot3(noise_freqs, noise_weights, ones(1,length(noise_freqs)).*max(abs(val(:)))*1.1 ,'k');
%  hold off;
% 
% 
% W = -mag2db(Quiet_Sample__Weight_Vs_Frequency);
% W_ = permute( Tools.interpVal_2D(W, Frequencies, Weights, noise_freqs, noise_weights, 'spline'), [2 1]);
% W_=W_(:);
% 
% %hold on
% figure(2);
% plot(noise_freqs, val/2);
% set(gca,'XScale','log');
% set(gca,'YScale','lin');
% set(gca,'TickDir','out');
% ylabel('Mag Difference');
% xlabel('Frequency (Hz)');
% %ylim([0 1]);
% grid on;
% 
% hold on
% %figure(3);
% plot(noise_freqs, W_);
% set(gca,'XScale','log');
% set(gca,'YScale','lin');
% set(gca,'TickDir','out');
% ylabel('Level (dB)');
% xlabel('Frequency (Hz)');
% grid on;