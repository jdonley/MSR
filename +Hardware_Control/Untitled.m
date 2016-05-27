[data,fs] = audioread('Z:\_RECOR~1\_PHYSI~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_1DB3B~1\Female_SA2_Original.WAV');
load('Z:\_RECOR~1\_PHYSI~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_1DB3B~1\Female_SA2_Quiet.mat')
load('Z:\_RECOR~1\_PHYSI~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_1DB3B~1\Female_SA2_Bright.mat')
orig = data;
Rec_Bright = double(Rec_Sigs_B');
Rec_Quiet = double(Rec_Sigs_Q');
Fs = 48000;
signal_info.Fs = 16000;
if Fs ~= signal_info.Fs
down_rate = Fs / signal_info.Fs ;
Rec_Bright_down = zeros(ceil(size(Rec_Bright).*[1 1/down_rate]));
Rec_Quiet_down = zeros(ceil(size(Rec_Quiet).*[1 1/down_rate]));
for  r = 1:size(Rec_Bright,1)
Rec_Bright_down(r,:) = ...
decimate( Rec_Bright(r,:), down_rate );
Rec_Quiet_down(r,:) = ...
decimate( Rec_Quiet(r,:), down_rate );
end
Rec_Bright = Rec_Bright_down;
Rec_Quiet  = Rec_Quiet_down ;
if ~isempty(orig)
orig = decimate( orig, down_rate );
end
end
orig(length(orig):size(Rec_Bright,2))=0; % Resize the original signal because the reverberant signal will be longer
if (length(orig) ~= length(Rec_Bright)) || (length(orig) ~= length(Rec_Quiet))
error('Size of the original signal does not match the reproduced signal!');
end
%c_speed = 343;%343m/s speed of sound in air
%max_delay = speaker_radius*2 / c_speed * signal_info.Fs;
max_delay = signal_info.Fs / 2;
Original = zeros(size(Rec_Bright,1),2,length(orig));
for r = 1:size(Rec_Bright,1)
delay = sigalign(Rec_Bright(r,:), orig, [-1 1]*max_delay) - 1;
if delay <= 0
Original(r,1,:) = [orig(-delay:end); zeros(-delay-1,1)];
elseif delay>0
Original(r,1,:) = orig;
Rec_Bright(r,:) = [Rec_Bright(r,delay:end), zeros(1,delay-1)];
end
delay = sigalign( Rec_Quiet(r,:), orig, [-1 1]*max_delay) - 1;
if delay <= 0
Original(r,2,:) = [orig(-delay:end); zeros(-delay-1,1)];
elseif delay>0
Original(r,2,:) = orig;
Rec_Quiet(r,:) = [Rec_Quiet(r,delay:end), zeros(1,delay-1)];
end
end
Tools.pesq_mex_vec(orig,Rec_Bright,16000,1)