% Microphone = dsp.AudioRecorder;
% Speaker = dsp.AudioPlayer;
% SpecAnalyzer = dsp.SpectrumAnalyzer;
% 
% %% Create filter
% fs = Microphone.SampleRate;
% [b,a] = cheby1(6,1,[1000]/(fs/2),'low');
% h = impz(b,a,Microphone.SamplesPerFrame);
% h_=repmat(h,1,2);
% %% Set hardware
% Microphone.DeviceName = Microphone.set.DeviceName{:}{3};
% Microphone.DeviceName = 'Stereo Mix (Realtek High Definition Audio)'
% Speaker.DeviceName = Speaker.set.DeviceName{:}{3};
% 
% %%
% while true
% audio = step(Microphone);
% % audio = Tools.fconv(audio,h_);    
% step(SpecAnalyzer,audio);
% % Speaker.step(audio./max(abs(audio(:))));
% end



%% Read from File and Write to Audio Device
% Read an MP3 audio file and play it through your default audio output
% device.

%%
% Create a |dsp.AudioFileReader| System object(TM) with default settings.
% Use the |audioinfo| function to return a structure containing information
% about the audio file.
fileReader = dsp.AudioFileReader('speech_dft.mp3','SamplesPerFrame',512);
fileInfo = audioinfo('speech_dft.mp3');

%%
% Create an |audioDeviceWriter| System object and specify the sample rate.
% Call |setup| to reduce the computational load of initialization in an
% audio stream loop.
deviceWriter = audioDeviceWriter(...
    'SampleRate',fileInfo.SampleRate, ...
    'SupportVariableSizeInput', true, ...
    'BufferSize', 1024);
setup(deviceWriter,...
    zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels));

%% Create filter
fs = fileReader.SampleRate;
[b,a] = cheby1(6,1,[1000]/(fs/2),'low');
% h = impz(b,a,fileReader.SamplesPerFrame);

%%
% In an audio stream loop, read an audio signal frame from the file, and
% write the frame to your device.
ad=[];
while ~isDone(fileReader)
    audioData = step(fileReader);
    audioData = filter(b,a,audioData);
    ad=[ad;audioData];   plot(ad);drawnow;
    play(deviceWriter,audioData);
end

%%
% Close the input file and release the device.
release(fileReader);
release(deviceWriter);