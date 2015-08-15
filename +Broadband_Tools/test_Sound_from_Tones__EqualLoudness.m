frequencies = linspace(100, 8000, 1024);
Fs = 16000;
levels = Perceptual_Tools.Loudness(frequencies, 60);
levels = db2mag(levels);
nSec = 100000/Fs;
y = 0;
for f = 1:length(frequencies)
    y = y + levels(f)*sin(linspace(0,nSec * frequencies(f) * 2 * pi + rand*2*pi, round(nSec*Fs)));
end
%%
%sound(x,Fs);
sound(y((end*3/4):end), Fs); 