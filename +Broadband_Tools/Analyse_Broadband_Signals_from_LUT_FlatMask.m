function Analyse_Broadband_Signals_from_LUT_FlatMask( Input_file_path, Input_file_name, Input_file_ext, LUT_resolution, Noise_Mask_dB, weight, loudspeaker_setup )
%Analyse_Broadband_Signals_from_LUT_FlatMask Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables
%Input_file_path = '+Miscellaneous\'; % Can be relative or exact
%Input_file_name = 'MaleSpeech16k';
%Input_file_ext  = '.wav';
% if nargin < 5
%     LUT_resolution = '256f_128w';
% elseif nargin < 6
%     LUT_resolution = [num2str(N_LUT_frequencies) 'f_128w'];
% elseif nargin < 7
% 	LUT_resolution = [num2str(N_LUT_frequencies) 'f_' num2str(N_LUT_weights) 'w'];
% end

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

loudspeakers   = loudspeaker_setup;  % Number of loudspeakers
if loudspeaker_setup == 16
    speaker_arc    = 180;  % Degrees
else
    speaker_arc    = 360;  % Degrees
end
first_speaker  = 0; % Degrees
speaker_radius = 1.5; % Metres

Output_file_path     = '+Speaker_Signals\'; % Can be relative or exact
Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
SetupInfo            = ['_' num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz_' ...
    num2str(angle_pw) 'pwAngle_' ...
    num2str(Noise_Mask_dB) 'dB_' ...
    num2str(weight) 'weight__withFlatMask'];
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
Nplanewaves = round( +Tools.interpVal(ceil((300-30)/((f_samples-1)^4)*(0:(f_samples-1)).^4+30), ...
    logspace(log10(150),log10(8000),f_samples), ...
    Frequencies_));

%Loudspeaker_Weights = zeros(f_samples, loudspeakers);
Bright_Zone_sample  = zeros(f_samples, 1);
Quiet_Zone_sample   = zeros(f_samples, 1);










%% Secondly, build a flat spectra desired multizone soundfield for all frequencies from the previous fft and save the speaker signals for each frequency bin.

load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat']);

LUT_Bright = Bright_Sample__Weight_Vs_Frequency;%(db2mag(Contrast__Weight_Vs_Frequency).^0.5 .* db2mag(Quiet_SPL__Weight_Vs_Frequency-94));
LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;%db2mag(Quiet_SPL__Weight_Vs_Frequency-94);

Bright_Zone_sample = permute( Tools.interpVal_2D(LUT_Bright, Frequencies, Weights, Frequencies_, weight), [2 1]);
Quiet_Zone_sample  = permute( Tools.interpVal_2D(LUT_Quiet, Frequencies, Weights, Frequencies_, weight), [2 1]);

if ~exist(['+Soundfield_Database\+From_LUT\' Output_file_path_ext],'dir'); mkdir(['+Soundfield_Database\+From_LUT\' Output_file_path_ext]); end
save(['+Soundfield_Database\+From_LUT\' Output_file_path_ext 'Weights_and_Samples__' SetupInfo '.mat'], ...
    ...%'Loudspeaker_Weights', ...
    'Bright_Zone_sample', ...
    'Quiet_Zone_sample' );










%% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal
load(['+Soundfield_Database\+From_LUT\' Output_file_path_ext 'Weights_and_Samples__' SetupInfo '.mat']);

% % Here we want to build the speaker signals for each speaker so that our loudspeaker weights are taken into account.
% % We want to form the entire spectrum by adding the conjugate of the frame
% % to the existing frame where the negative frequencies of the transform
% % would usually exist.
% Loudspeaker_Weights = [zeros(trunc_index_low-1, loudspeakers); Loudspeaker_Weights; zeros( size(Z,2) - trunc_index_high, loudspeakers) ];
% Loudspeaker_Weights = permute( repmat(Loudspeaker_Weights, [1 1 size(Z,1)]), [3 1 2]);
% Z = repmat(Z, [1 1 loudspeakers]);
% 
Loudspeakers_ = zeros( [size(Z,1) (size(Z,2))*2 size(Z,3)] );
% for spkr = 1:loudspeakers
%     Loudspeakers_(:,1:end/2,spkr) = Z(:,:,spkr) .* Loudspeaker_Weights(:,:,spkr);
%     Loudspeakers_(:,end/2+1:end,spkr) = conj( [-Loudspeakers_(:,1,spkr).*0 Loudspeakers_(:,end/2:-1:2,spkr)]);
% end
Original = [Z(:,:,1) conj( [-Z(:,1,1).*0 Z(:,end:-1:2,1)] )];
% clear Loudspeaker_Weights; % Save on memory
% 
% % We then want to perform an Inverse FFT (ifft) on each full spectrum frame
 for frame = 1:size(Loudspeakers_, 1)
%     for spkr = 1:loudspeakers
%         Loudspeakers_(frame,:,spkr) = ifft( Loudspeakers_(frame,:,spkr) );
%     end
     Original(frame,:) = ifft( Original(frame, :) );
 end
% 
% % Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
% %Loudspeaker_Signals = zeros( [(size(Z,1)+ceil(overlap))*size(Z,2)*2*(1-overlap) loudspeakers] );
% for spkr = 1:loudspeakers
%     Loudspeaker_Signals(:,spkr) = Tools.OverlapAdd( Loudspeakers_(:,:,spkr), overlap );
% end
 Original_ = Tools.OverlapAdd( Original, overlap );
% clear Loudspeakers_; % Save on memory
% 
% % Normalise Loudspeaker Signals
% Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)));
 Original_ = Original_ ./ max(abs(Original_(:)));
% 
% 
% % Once we have the speaker signals we should save them for later use as .wav files
% filenumbers = num2str((1:loudspeakers)');
% filenumbers(filenumbers==' ') = '_';
% fullpath = [repmat([Output_file_path Output_file_path_ext Output_file_name], [loudspeakers 1]) ...
%     filenumbers ...
%     repmat(Output_file_ext, [loudspeakers 1]) ];
% if ~exist([Output_file_path Output_file_path_ext],'dir'); mkdir([Output_file_path Output_file_path_ext]); end
% 
% for spkr = 1:loudspeakers
%     audiowrite(fullpath(spkr,:), Loudspeaker_Signals(:, spkr), Fs);
% end
% audiowrite([Output_file_path Output_file_path_ext ...
%     Input_file_name '__Original_' ...
%     num2str(f_low ) 'Hz-' ...
%     num2str(f_high) 'Hz' ...
%     Input_file_ext], Original_, Fs);










%%
Bright_Zone_sample = [zeros(trunc_index_low-1, 1); Bright_Zone_sample; zeros( size(Z,2) - trunc_index_high, 1) ];
Bright_Zone_sample = permute( repmat(Bright_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);

Quiet_Zone_sample  = [zeros(trunc_index_low-1, 1); Quiet_Zone_sample; zeros( size(Z,2) - trunc_index_high, 1) ];
Quiet_Zone_sample  = permute( repmat(Quiet_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);

Bright = zeros( [size(Z,1) (size(Z,2))*2] );
Quiet  = zeros( [size(Z,1) (size(Z,2))*2] );

Bright(:,1:end/2)     = squeeze( Z(:,:,1) ) .* Bright_Zone_sample(:,:);
Bright(:,end/2+1:end) = conj( [-Bright(:,1).*0 Bright(:,end/2:-1:2)] );
Quiet(:,1:end/2)      = squeeze( Z(:,:,1) ) .* Quiet_Zone_sample;
Quiet(:,end/2+1:end)  = conj( [-Quiet(:,1).*0 Quiet(:,end/2:-1:2)] );

% We then want to perform an Inverse FFT (ifft) on each full spectrum frame
for frame = 1:size(Z, 1)
    Bright(frame,:) = ifft( Bright(frame, :) );
    Quiet( frame,:) = ifft( Quiet( frame, :) );
end

% Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
Bright_Signal = Tools.OverlapAdd( Bright, overlap );
Quiet_Signal  = Tools.OverlapAdd( Quiet , overlap );

%Add flat noise mask
Masking_Noise = max(Bright_Signal)/2 * (rand(size(Bright_Signal))*2 -1) * db2mag( - Noise_Mask_dB ); % Add (- Noise_Mask_dB) dB of noise to the signal
Bright_Signal = Bright_Signal + Masking_Noise;
Quiet_Signal  = Quiet_Signal  + Masking_Noise;

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

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_Results_' num2str(weight) 'weight_withFlatMask.csv'],'a');
PESQ_Bright = pesq(Original_path,Bright_path);
PESQ_Quiet = pesq(Original_path,Quiet_path);
fprintf(fileID,'%s,PESQ_Bright,%f,PESQ_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_Bright, PESQ_Quiet);
fclose(fileID);

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_MOS_Results_' num2str(weight) 'weight_withFlatMask.csv'],'a');
PESQ_MOS_Bright = pesq2mos(PESQ_Bright);
PESQ_MOS_Quiet = pesq2mos(PESQ_Quiet);
fprintf(fileID,'%s,PESQ_MOS_Bright,%f,PESQ_MOS_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_MOS_Bright, PESQ_MOS_Quiet);
fclose(fileID);

%% Calculate and save Active Speech Level values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'ASL_Results_' num2str(weight) 'weight_withFlatMask.csv'],'a');

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

fileID = fopen(['+Results\' Output_file_path_ext 'SI_Results_' num2str(weight) 'weight_withFlatMask.csv'],'a');

d_SI_Bright = Tools.stoi(Original_, Bright_Signal, Fs);
d_SI_Quiet  = Tools.stoi(Original_, Quiet_Signal,  Fs);

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


WordsCorrect_SI_Bright = Tools.stoi_d2percCorr(d_SI_Bright);
WordsCorrect_SI_Quiet  = Tools.stoi_d2percCorr(d_SI_Quiet);

fprintf(fileID,'%s,d_SI_Bright,%f,d_SI_Quiet,%f,WordsCorrect_SI_Bright,%f,WordsCorrect_SI_Quiet,%f,\r\n', ...
        [Input_file_name SetupInfo], ...
        d_SI_Bright, ...
        d_SI_Quiet, ...
        WordsCorrect_SI_Bright, ...
        WordsCorrect_SI_Quiet );
fclose(fileID);

end