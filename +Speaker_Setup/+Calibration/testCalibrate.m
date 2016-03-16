clc;clear;
delete(gcp);
current_pool = parpool; %Start new pool

%%
Nspkrs = 24;
Dur = 40; %seconds
Fs = 48000; %Hz
PadDur = [5, 5]; %seconds
F_band = [20, 20000]; %Hz
IR_len = 0.5; %seconds
Reg = [20 -6]; %[passband_gain, stopband_gain] (dB)

%%
fprintf('\nReading audio files... \t\t\t\t\t');
y = audioread('Z:\+Recordings\Chamber\Rec2_4-3-16_11-30pm.wav');
y_in = audioread('Z:\+Recordings\Chamber\sinesweep_20Hz-20kHz_30sec_5secTail_5secHead.WAV');


fprintf('Done\n');
%%
fprintf('Segmenting loudspeaker recordings ... \t');
[R, y_inNoPad] = Speaker_Setup.Calibration.SplitRecording( y, y_in, Nspkrs, Dur, Fs );


fprintf('Done\n');
%%
fprintf('Determining impulse responses ... \t\t');
IRs = Speaker_Setup.Calibration.IRsFromRecordings( R, y_inNoPad, F_band, Fs );


fprintf('Done\n');
%%
fprintf('Computing inverse impulse responses ... ');
invIRs = Speaker_Setup.Calibration.getInverseFilters( IRs, IR_len, F_band, Fs, Reg, 'toeplitz' );


fprintf('Done\n');
%%
fprintf('->\tEqualising frequency response ... \t');
irs = Tools.fconv( IRs, invIRs );
 

fprintf('Done\n');
%%
fprintf('Computing gain adjustment ... \t\t\t');
IRgains = Speaker_Setup.Calibration.getGainAdjustments( irs, F_band, Fs );


fprintf('Done\n');
%%
fprintf('Computing time-alignment ... \t\t\t');
TimeAlignIRs = Speaker_Setup.Calibration.getTimeAlignIRs( IRs );


fprintf('Done\n');
%%
fprintf('Creating filters ... \t\t\t\t\t');
filts = Tools.fconv( ...
    invIRs .* repmat( IRgains, [size(invIRs,1), 1]), ...
    TimeAlignIRs );


fprintf('Done\n');
%%
fprintf('\nFinished\n');