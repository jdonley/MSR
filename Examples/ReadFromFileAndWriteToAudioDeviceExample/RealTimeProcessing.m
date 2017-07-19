frameLength = 256;
fileReader = dsp.AudioFileReader(...
    'Counting-16-44p1-mono-15secs.wav',...
    'SamplesPerFrame',frameLength);
deviceWriter = audioDeviceWriter(...
    'SampleRate',fileReader.SampleRate);
scope = dsp.TimeScope(...
    'SampleRate',fileReader.SampleRate,...
    'TimeSpan',16,...
    'BufferLength',1.5e6,...
    'YLimits',[-1,1]);
dRG = noiseGate(...                                         %<---
    'SampleRate',fileReader.SampleRate,...                  %<---
    'Threshold',-25,...                                     %<---
    'AttackTime',10e-3,...                                  %<---
    'ReleaseTime',20e-3,...                                 %<---
    'HoldTime',0);                                          %<---

while ~isDone(fileReader)
    signal = fileReader();
    noisySignal = signal + 0.0025*randn(frameLength,1);     %<---
    processedSignal = dRG(noisySignal);                     %<---
    deviceWriter(processedSignal);                          %<---
    scope([noisySignal,processedSignal]);                   %<---
end

release(fileReader);
release(deviceWriter);
release(scope);
release(dRG);                                               %<---