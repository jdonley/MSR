clc;clear;

%%
frequency_range = [50, 20000]; %Hz
regularisation = [60 -6]; %[passband_gain, stopband_gain] (dB)
filter_length = 0.5; %seconds

[reference, fs] = audioread('Z:\+Recordings\Chamber\sinesweep_20Hz-20kHz_30sec_5secTail_5secHead.WAV');

%% Calibrated
% recording = audioread('Z:\+Recordings\Chamber\Rec1_8-3-16_12-10am__calibrated.wav');
% dur = length(audioread('Z:\+Speaker_Signals\CalibratedSweeps\Loudspeaker_1.wav'))/fs; %seconds

%% Non-calibrated
recording = audioread('Z:\+Recordings\Chamber\Rec2_4-3-16_11-30pm.wav');
dur = length(reference) / fs; %seconds

%%
Nspkrs = floor(length(recording) / dur / fs );

%Segment loudspeaker recordings
[R, referenceSigNoPad] = Speaker_Setup.Calibration.SplitRecording( recording, reference, Nspkrs, dur, fs );

%Determine impulse responses
IRs = Speaker_Setup.Calibration.IRsFromRecordings( R, referenceSigNoPad, frequency_range, fs );


%% Plot
figure(1);
periodogram(IRs,[],1024,fs);
set(gca,'XScale','log');
grid minor;
axis([0.1, 10, -180, -140]);

%%
figure(2);
timeMat = repmat( linspace(0,size(IRs,1)/fs   * 1000,size(IRs,1))', [1 Nspkrs]);

[pkval,pk]=max(IRs);
centPk = max(pk-1)/fs  * 1000;

plot((timeMat - centPk)*343,IRs);
grid minor;
hold on;
plot(((pk-1)/fs  * 1000 - centPk)* 343,pkval,'or');
hold off;

xlim([-1 1] .* 100);
xlabel('Millimeters')
ylabel('Amplitude')
title('Impulse Response from 24 Loudspeakers');
text(20,10e-4,'Sampling Frequency: 48kHz');