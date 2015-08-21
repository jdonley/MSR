function Broadband_PrivacyWeighted_v5_Mixed_ANC_Mask( Input_file_path, Input_file_name, Input_file_ext, LUT_resolution, Noise_Mask_dB, loudspeaker_setup)
%Analyse_Broadband_Signals_from_LUT Summary of this function goes here
%   Detailed explanation goes here
Version = 'v5';

%% Setup Variables
%Input_file_path = '+Miscellaneous\'; % Can be relative or exact
%Input_file_name = 'MaleSpeech16k';
%Input_file_ext  = '.wav';
%LUT_resolution = '256f_128w';

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
quiet_mask_noise_dB = Noise_Mask_dB; % Phons (dB)

if loudspeaker_setup == 16
    loudspeakers   = 16;  % Number of loudspeakers
    speaker_arc    = 180;  % Degrees
elseif loudspeaker_setup == 65
    loudspeakers   = 65;  % Number of loudspeakers
    speaker_arc    = 360;  % Degrees
else
    error('Number of loudspeakers not currently supported');
end
first_speaker  = 90+15; % Degrees
speaker_radius = 1.5; % Metres


Output_file_path     = '+Speaker_Signals\'; % Can be relative or exact
Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
SetupInfo            = ['_' num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz_' ...
    num2str(angle_pw) 'pwAngle_' ...
    num2str(Noise_Mask_dB) 'dB_privacy_' Version '__'];
Output_file_name     = [Input_file_name '__' ...
    num2str(loudspeakers) 'spkrs_' ...
    SetupInfo];
Output_file_ext      = '.wav';





Bright_Signal = [];
Quiet_Signal  = [];
ANC_Bright_Signal = [];
ANC_Quiet_Signal  = [];

%% For ANC we want to produce the resulting signal of the quiet zone back into the quiet zone but inverted

[Input_Signal, Fs_file] = audioread( [Input_file_path Input_file_name Input_file_ext] );

for ANC_Step = 1:2
    
    %% Firstly, find the frequency domain representation of an audio file that is wished to be reproduced in the spatial domain.
    [Z, Frequencies_, ~, Window] = Broadband_Tools.FFT_custom_vec(Input_Signal, Nfft, Fs, overlap);
    
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
    %For the first step of ANC we have the zones set to how we want them
    %For the second step we switch the zones so we can generate the cancelling signal    
    if ANC_Step == 1
        load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat']);
    elseif ANC_Step == 2
        load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '_zones_swapped.mat']);
    end    
    
    LUT        = Quiet_SPL__Weight_Vs_Frequency(: , 1:end) - 94;%... %For each frequency set
    LUT_Bright = Bright_Sample__Weight_Vs_Frequency;
    LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;
    
    %We need to find a set of weights for each frame that makes each individual SPL suppress to a
    %given value (This values is: quiet_mask_noise_dB)
    SPL_values = 20 * log10( abs( Z(:,trunc_index_low:trunc_index_high) ) .^ (1/2) * 10^(SPL_reference / 20));
    
    % Find the weights that provide attenuation just below the threshold
    %First find the desired SPLs
    %Desired_SPL = Perceptual_Tools.Threshold_in_Quiet( Frequencies_in_zones(1, :) , 'ISO226'); % ISO226 Standard
    Desired_SPL = Perceptual_Tools.Loudness( Frequencies_ , 40);
    Desired_SPL = repmat(Desired_SPL, [1 size(SPL_values,1)])';
    
    % Because the frequencies in the LUT don't match up we will interpolate
    % We will linearly interpolate the closest index that will gives us a weight that
    % provides attenuation at the threshold
    SPL_Difference = SPL_values - Desired_SPL;
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
    
    
    
    
    %% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal
    
    Loudspeakers_ = zeros( [size(Z,1) (size(Z,2))*2 size(Z,3)] );
    
    Original = [Z(:,:,1) conj( [-Z(:,1,1).*0 Z(:,end:-1:2,1)] )];
    
    for frame = 1:size(Loudspeakers_, 1)
        Original(frame,:) = ifft( Original(frame, :) );
    end
    
    Original_ = +Broadband_Tools.OverlapAdd( Original, overlap );
    
    Original_ = Original_ ./ max(abs(Original_(:)));
    
    if ANC_Step == 1
        Original_Input_Signal = Original_;
    end
    
    
    
    
    
    
    
    %%
    Bright_Zone_sample = [zeros( size(Z,1), trunc_index_low-1) Bright_Zone_sample zeros( size(Z,1), size(Z,2) - trunc_index_high) ];
    %Bright_Zone_sample = permute( repmat(Bright_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);
    
    Quiet_Zone_sample  = [zeros( size(Z,1), trunc_index_low-1) Quiet_Zone_sample zeros( size(Z,1), size(Z,2) - trunc_index_high) ];
    %Quiet_Zone_sample  = permute( repmat(Quiet_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);
    
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
%     Bright_sig = overlapadd(Bright, Window, length(Window) * overlap);
%     Quiet_sig = overlapadd(Quiet, Window, length(Window) * overlap);
    Bright_sig = Broadband_Tools.OverlapAdd( Bright, overlap );
    Quiet_sig  = Broadband_Tools.OverlapAdd( Quiet , overlap );
    
    % Normalise Signals
    maxVal = max(abs(Bright_sig(:)));
    Bright_sig = Bright_sig ./ maxVal;
    Quiet_sig  = Quiet_sig  ./ maxVal;
    
    if ANC_Step == 1
        Bright_Signal = Bright_sig;
        Quiet_Signal = Quiet_sig;
    elseif ANC_Step == 2
        ANC_Bright_Signal = -Bright_sig;
        ANC_Quiet_Signal = -Quiet_sig;
    end

    % For ANC we need to flip the zones and make our next input signal the
    % inverted quiet zone result which should cancel the previous quiet
    % zone (spatial errors will cause this to not work perfectly)
    % This may also cause some unwanted cancellation in the bright zone,
    % however it should be suppressed
    Masking_Noise = max(Quiet_sig)/2 * (rand(size(Quiet_sig))*2 -1) * db2mag( - Noise_Mask_dB ); % Add (- Noise_Mask_dB) dB of noise to the signal
    Input_Signal = Quiet_sig + Masking_Noise; % Parallel Active Noise Cancellation
    
end




%% We then need to normalise the ANC Quiet Signal to the same level as the previous quiet signal
% To do this we need to find the cross-correlation peak and divide it by
% the energy in the first Quiet signal to find the scaling factor
% The signals are assumed here to be alligned in the time-domain

Quiet_corr = xcorr( Quiet_Signal, ANC_Bright_Signal);
peakVal = max(abs(Quiet_corr));
Energy = sum( Quiet_Signal .^ 2 );

if isnan(peakVal)
    peakVal = mean([max(abs(Quiet_Signal)), max(abs(ANC_Bright_Signal))]);
end

Scaling_Factor = peakVal / Energy;

ANC_Bright_Signal = ANC_Bright_Signal * 1/Scaling_Factor;
ANC_Quiet_Signal  = ANC_Quiet_Signal * 1/Scaling_Factor;

Bright_Signal = Bright_Signal + ANC_Quiet_Signal;
Quiet_Signal  = Quiet_Signal  + ANC_Bright_Signal;

% Normalise Loudspeaker Signals
maxVal = max(abs(Bright_Signal(:)));
Bright_Signal = Bright_Signal ./ maxVal;
Quiet_Signal  = Quiet_Signal  ./ maxVal;



%%

% Once we have the speaker signals we should save them for later use as .wav files
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end
Original_path = ['+Results\' Output_file_path_ext ...
    Input_file_name '__Original_' ...
    num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz' ...
    Input_file_ext];
audiowrite(Original_path, Original_Input_Signal, Fs);
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

% Wait for files to write because they are used soon after
% Tools.wait_for_file(Original_path);
% Tools.wait_for_file(Bright_path);
% Tools.wait_for_file(Quiet_path);






%% Calculate and save PESQ Values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_Results_PrivacyWeighted_' Version '.csv'],'a');

PESQ_Bright = pesq(Original_path,Bright_path);
PESQ_Quiet = pesq(Original_path,Quiet_path);

fprintf(fileID,'%s,PESQ_Bright,%f,PESQ_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_Bright, PESQ_Quiet);
fclose(fileID);

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_MOS_Results_PrivacyWeighted_' Version '.csv'],'a');
PESQ_MOS_Bright = pesq2mos(PESQ_Bright);
PESQ_MOS_Quiet = pesq2mos(PESQ_Quiet);
fprintf(fileID,'%s,PESQ_MOS_Bright,%f,PESQ_MOS_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_MOS_Bright, PESQ_MOS_Quiet);
fclose(fileID);


%% Calculate and save Active Speech Level values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'ASL_Results_PrivacyWeighted_' Version '.csv'],'a');

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

fileID = fopen(['+Results\' Output_file_path_ext 'SI_Results_PrivacyWeighted_' Version '.csv'],'a');

d_SI_Bright = Tools.stoi(Original_Input_Signal, Bright_Signal, Fs);
d_SI_Quiet  = Tools.stoi(Original_Input_Signal, Quiet_Signal,  Fs);

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