function Broadband_PrivacyWeighted( Input_file_path, Input_file_name, Input_file_ext, Noise_Mask_dB)
%Analyse_Broadband_Signals_from_LUT Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables
%Input_file_path = '+Miscellaneous\'; % Can be relative or exact
%Input_file_name = 'MaleSpeech16k';
%Input_file_ext  = '.wav';
LUT_resolution = '512f_256w';

Fs = 16000; % Sampling frequency
Nfft = 1024;% Number of fft components
overlap = 0.5;
f_low  = 150;  % Hz
f_high = 8000; % Hz

resolution = 100;
phase = 0;
radius = 0.3;
angle_pw = 15;
radius2origin = 0.6;
angle2origin  = [0 180];
%weight = 0.05;
Weights_ = [];

SPL_reference = 94; % 94dB is 1Pa (RMS) when reference to 20uPa (RMS)
quiet_mask_noise_dB = 30;%Noise_Mask_dB; % Phons (dB)

loudspeakers   = 65;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
first_speaker  = 90+15; % Degrees
speaker_radius = 1.5; % Metres

Output_file_path     = '+Speaker_Signals\'; % Can be relative or exact
Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
SetupInfo            = [num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz_' ...
    num2str(angle_pw) 'pwAngle_' ...
    num2str(Noise_Mask_dB) 'dB_privacy__'];
Output_file_name     = [Input_file_name '__' ...
    num2str(loudspeakers) 'spkrs_' ...
    SetupInfo];
Output_file_ext      = '.wav';










%% Firstly, find the frequency domain representation of an audio file that is wished to be reproduced in the spatial domain.
[Z, Frequencies_] = Broadband_Tools.FFT_custom([Input_file_path Input_file_name Input_file_ext], Nfft, Fs, overlap);

% Truncate to frequencies in the range f_low <-> f_high
trunc_index_low  = find(Frequencies_ < f_low , 1, 'last' ) + 1;
trunc_index_high = find(Frequencies_ > f_high, 1 ) + 1;
if isempty(trunc_index_low)
    trunc_index_low = 1;
end
if isempty(trunc_index_high)
    trunc_index_high = length(Frequencies_);
end
%Z = Z ( :, trunc_index_low:trunc_index_high );
Frequencies_ = Frequencies_( :, trunc_index_low:trunc_index_high );

% Work out the coresponding number of planewaves for these frequencies
f_samples = length(Frequencies_);
frames    = size(Z,1);
Nplanewaves = round( +Tools.interpVal(ceil((300-30)/((f_samples-1)^4)*(0:(f_samples-1)).^4+30), ...
    logspace(log10(150),log10(8000),f_samples), ...
    Frequencies_));

%Loudspeaker_Weights = zeros(f_samples, loudspeakers);
Bright_Zone_sample  = zeros(f_samples, 1);
Quiet_Zone_sample   = zeros(f_samples, 1);










%% Secondly, build a flat weighted desired multizone soundfield for all frequencies from the previous fft and save the bright and quiet samples for each frequency bin.

load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat']);

LUT        = Quiet_SPL__Weight_Vs_Frequency(: , 1:end) - 94;%... %For each frequency set
LUT_Bright = Bright_Sample__Weight_Vs_Frequency;
LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;

%We need to find a set of weights for each frame that makes each individual SPL suppress to a
%given value (This values is: quiet_mask_noise_dB)
SPL_values = 20 * log10( abs( Z(:,trunc_index_low:trunc_index_high) ) .^ (1/2) * 10^(SPL_reference / 20));

% Find the weights that provide attenuation just below the threshold
%First find the desired SPLs
%Desired_SPL = Perceptual_Tools.Threshold_in_Quiet( Frequencies_in_zones(1, :) , 'ISO226'); % ISO226 Standard
Desired_SPL = Perceptual_Tools.Loudness( Frequencies_ , quiet_mask_noise_dB);
Desired_SPL = repmat(Desired_SPL, [1 size(SPL_values,1)])';
% Alternative Desired_SPL
% lut_max = max(LUT(:,trunc_index_low:trunc_index_high),[],1);
% lut_min = min(LUT(:,trunc_index_low:trunc_index_high),[],1);
% lut_control = (lut_max - lut_min) * 4/4;
% Desired_SPL = imfilter(SPL_values,fspecial('average', 20 ),'replicate') - repmat(lut_control,[size(SPL_values,1) 1]);
% Desired_SPL = Desired_SPL + (rand(size(Desired_SPL))-0.5) * 30;

% Because the frequencies in the LUT don't match up we will interpolate
% We will linearly interpolate the closest index that will gives us a weight that
% provides attenuation at the threshold
SPL_Difference = Desired_SPL - SPL_values ;

Weights_ = zeros([f_samples frames]);
for frame = 1:frames    
Weights_(:,frame) = Tools.interpFromVal_2D(LUT, Frequencies, Weights, Frequencies_, SPL_Difference(frame,:));
end
%surf(1:frames,Frequencies_,Weights_,'LineStyle','none');set(gca, 'ZScale', 'log');set(gca, 'YScale', 'log');view(2);

% We then choose the coressponding value in the LUT and gather the samples
% associated with it
Bright_Zone_sample = zeros([frames f_samples]);
Quiet_Zone_sample  = zeros([frames f_samples]);
for frame = 1:frames
Bright_Zone_sample(frame, :) = permute( Tools.interpVal_2D(LUT_Bright, Frequencies, Weights, Frequencies_, Weights_(:,frame)'), [2 1]);
Quiet_Zone_sample(frame, :)  = permute( Tools.interpVal_2D(LUT_Quiet , Frequencies, Weights, Frequencies_, Weights_(:,frame)'), [2 1]);
end

% if quiet_mask_noise_dB == 40
%     figure(1);
%     plot(SPL_values(35,:)); hold on; plot(Desired_SPL(35,:),'-k');plot(SPL_Difference(35,:),'-r');plot(mag2db(abs(Bright_Zone_sample(35,:))) - mag2db(abs(Quiet_Zone_sample(35,:))) + Desired_SPL(35,:),'-g');set(gca, 'XScale', 'log');
%     figure(2);
%     hold on; plot(mag2db(abs(Bright_Zone_sample(35,:))),'-g');plot(mag2db(abs(Quiet_Zone_sample(35,:))),'-r');set(gca, 'XScale', 'log');
% end


%% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal

Loudspeakers_ = zeros( [size(Z,1) (size(Z,2))*2 size(Z,3)] );

Original = [Z(:,:,1) conj( [-Z(:,1,1).*0 Z(:,end:-1:2,1)] )];

 for frame = 1:size(Loudspeakers_, 1)
     Original(frame,:) = ifft( Original(frame, :) );
 end
 
 Original_ = +Broadband_Tools.OverlapAdd( Original, overlap );
 
 Original_ = Original_ ./ max(abs(Original_(:)));










%%
Bright_Zone_sample = [zeros( size(Z,1), trunc_index_low-1) Bright_Zone_sample zeros( size(Z,1), size(Z,2) - trunc_index_high) ];
%Bright_Zone_sample = permute( repmat(Bright_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);

Quiet_Zone_sample  = [zeros( size(Z,1), trunc_index_low-1) Quiet_Zone_sample zeros( size(Z,1), size(Z,2) - trunc_index_high) ];
%Quiet_Zone_sample  = permute( repmat(Quiet_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);

Bright = zeros( [size(Z,1) (size(Z,2))*2] );
Quiet  = zeros( [size(Z,1) (size(Z,2))*2] );

Bright(:,1:end/2)     = squeeze( Z(:,:,1) ) .* Bright_Zone_sample;%abs(Bright_Zone_sample) .* exp( 1i * angle(squeeze( Z(:,:,1) )) );
Bright(:,end/2+1:end) = conj( [-Bright(:,1).*0 Bright(:,end/2:-1:2)] );
Quiet(:,1:end/2)      = squeeze( Z(:,:,1) ) .* Quiet_Zone_sample;%abs(Quiet_Zone_sample) .* exp( 1i * angle(squeeze( Z(:,:,1) )) );
Quiet(:,end/2+1:end)  = conj( [-Quiet(:,1).*0 Quiet(:,end/2:-1:2)] );

% if quiet_mask_noise_dB == 40
%     figure(3); hold on; plot(mag2db(abs(Quiet(35,1:end/2))),'-g');plot(mag2db(abs(squeeze( Z(35,:,1) ))),'-r');plot(mag2db(abs(Quiet_Zone_sample(35,:))));set(gca, 'XScale', 'log');
% end


% plot(Frequencies_,Desired_SPL(35,:)); hold on;
% plot(Frequencies_,SPL_values(35,:));
% plot(Frequencies_,SPL_Difference(35,:));
% plot(Frequencies_,mag2db(abs(Quiet_Zone_sample(35,trunc_index_low:trunc_index_high))));
% plot(Frequencies_,mag2db(abs(Quiet(35,trunc_index_low:trunc_index_high))));  hold off;
% legend 'Desired SPL' 'SPL values' 'SPL Difference' 'Quiet Zone sample' 'Quiet'

%%
% We then want to perform an Inverse FFT (ifft) on each full spectrum frame
for frame = 1:size(Z, 1)
    Bright(frame,:) = ifft( Bright(frame, :) );
    Quiet( frame,:) = ifft( Quiet( frame, :) );
end

% Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
Bright_Signal = +Broadband_Tools.OverlapAdd( Bright, overlap );
Quiet_Signal  = +Broadband_Tools.OverlapAdd( Quiet , overlap );

% Normalise Loudspeaker Signals
maxVal = max(abs(Bright_Signal(:)));
Bright_Signal = Bright_Signal ./ maxVal;
Quiet_Signal  = Quiet_Signal  ./ maxVal;


% Once we have the speaker signals we should save them for later use as .wav files
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end
Original_path = ['+Results\' Output_file_path_ext ...
    Input_file_name '__Original_' ...
    num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz' ...
    Input_file_ext];
audiowrite(Original_path, Original_, Fs);
Bright_path = ['+Results\' Output_file_path_ext ...
    Input_file_name '___Bright_' ...
    SetupInfo ...
    Input_file_ext];
audiowrite(Bright_path, Bright_Signal, Fs);
Quiet_path = ['+Results\' Output_file_path_ext ...
    Input_file_name '___Quiet__' ...
    SetupInfo ...
    Input_file_ext];
audiowrite(Quiet_path, Quiet_Signal, Fs);










%% Calculate and save PESQ Values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_Results_Privacy_Weighted.csv'],'a');
PESQ_Bright = pesq(Original_path,Bright_path);
PESQ_Quiet = pesq(Original_path,Quiet_path);
fprintf(fileID,'%s,PESQ_Bright,%f,PESQ_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_Bright, PESQ_Quiet);
fclose(fileID);

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_MOS_Results_Privacy_Weighted.csv'],'a');
PESQ_MOS_Bright = pesq2mos(PESQ_Bright);
PESQ_MOS_Quiet = pesq2mos(PESQ_Quiet);
fprintf(fileID,'%s,PESQ_MOS_Bright,%f,PESQ_MOS_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_MOS_Bright, PESQ_MOS_Quiet);
fclose(fileID);


%% Calculate and save Active Speech Level values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'ASL_Results_Privacy_Weighted.csv'],'a');

[ASL_Bright_dB, Activity_Bright] = activlev(Bright_Signal, Fs, 'd');
[ASL_Quiet_dB,  Activity_Quiet ] = activlev(Quiet_Signal,  Fs, 'd');

fprintf(fileID,'%s,ASL_Bright_dB,%f,ASL_Quiet_dB,%f,Activity_Bright,%f,Activity_Quiet,%f,\r\n', ...
        [Input_file_name SetupInfo], ...
        ASL_Bright_dB, ...
        ASL_Quiet_dB, ...
        Activity_Bright, ...
        Activity_Quiet);
fclose(fileID);

%% Calculate and save Speech Intelligibility values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'SI_Results_Privacy_Weighted.csv'],'a');

d_SI_Bright = taal2011(Original_, Bright_Signal, Fs);
d_SI_Quiet  = taal2011(Original_, Quiet_Signal,  Fs);

% Following the publications by tall2011
%    References:
%      C. H. Taal, R. C. Hendriks, R. Heusdens, and J. Jensen. A Short-Time
%      Objective Intelligibility Measure for Time-Frequency Weighted Noisy
%      Speech. In Acoustics Speech and Signal Processing (ICASSP), pages
%      4214-4217. IEEE, 2010.
%      
%      C. H. Taal, R. C. Hendriks, R. Heusdens, and J. Jensen. An Algorithm
%      for Intelligibility Prediction of Time-Frequency Weighted Noisy Speech.
%      IEEE Transactions on Audio, Speech and Language Processing,
%      19(7):2125-2136, 2011.

% for IEEE English library
f_IEEE_a = -17.4906;
f_IEEE_b = 9.6921;

WordsCorrect_SI_Bright = 100 / (1 + exp(f_IEEE_a * d_SI_Bright + f_IEEE_b ) );
WordsCorrect_SI_Quiet  = 100 / (1 + exp(f_IEEE_a * d_SI_Quiet  + f_IEEE_b ) );

fprintf(fileID,'%s,d_SI_Bright,%f,d_SI_Quiet,%f,WordsCorrect_SI_Bright,%f,WordsCorrect_SI_Quiet,%f,\r\n', ...
        [Input_file_name SetupInfo], ...
        d_SI_Bright, ...
        d_SI_Quiet, ...
        WordsCorrect_SI_Bright, ...
        WordsCorrect_SI_Quiet );
fclose(fileID);

end