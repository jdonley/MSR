[x, Fs] = wavread('+Miscellaneous\crescendo.wav');
spec = abs(spectrogram(x));
[L_,f_] = findpeaks(spec(:,end/2),'NPeaks',1024,'MinPeakHeight',1,'MinPeakDistance',1);
stem(f_,L_);
nSec = length(x)/Fs;
y = 0;
for f = 1:length(f_)    
    y = y + L_(f)*sin(linspace(0,nSec * f_(f) * 2 * pi, round(nSec*Fs)));
end
%%
sound(x,Fs);
sound(y.*linspace(0.000001,0.0015,round(nSec*Fs)), Fs); 