%% Read from File and Write to Audio Device
% Read an MP3 audio file and play it through your default audio output
% device.

%%
% Create a |dsp.AudioFileReader| System object(TM) with default settings.
% Use the |audioinfo| function to return a structure containing information
% about the audio file.
fileReader = dsp.AudioFileReader('speech_dft.mp3');
fileInfo = audioinfo('speech_dft.mp3');

%%
% Create an |audioDeviceWriter| System object and specify the sample rate.
% Call |setup| to reduce the computational load of initialization in an
% audio stream loop.
deviceWriter = audioDeviceWriter(...
    'SampleRate',fileInfo.SampleRate);
setup(deviceWriter,...
    zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels));

%%
% In an audio stream loop, read an audio signal frame from the file, and
% write the frame to your device.
while ~isDone(fileReader)
    audioData = step(fileReader);
    play(deviceWriter,audioData);
end

%%
% Close the input file and release the device.
release(fileReader);
release(deviceWriter);