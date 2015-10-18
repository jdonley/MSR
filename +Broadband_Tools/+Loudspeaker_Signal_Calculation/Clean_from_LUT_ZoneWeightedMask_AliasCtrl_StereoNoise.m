function Clean_from_LUT_ZoneWeightedMask_AliasCtrl_StereoNoise( Input_file_path, Input_file_name, Input_file_ext, LUT_resolution, Noise_Mask_dB, weight, loudspeaker_setup, angle_pw, leakage_angle_delta )
%Clean_from_LUT_ZoneWeightedMask Summary of this function goes here

masker_type = 'ZoneWeightMaskAliasCtrlStereoNoise';
%% Setup Variables

Fs = 16000; % Sampling frequency
Nfft = 1024;% Number of fft components
overlap = 0.5;
f_low  = 150;  % Hz
f_high = 8000; % Hz

resolution = 100;
phase = 0;
radius = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leakage_angle = 0;
% STEREO NOISE
if nargin < 9
    leakage_angle_delta = 30; % Angle to virtual stereo noise source from angle of intereference (leakage)
end
leakage_angle_L = leakage_angle + leakage_angle_delta;
leakage_angle_R = leakage_angle - leakage_angle_delta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 8
    angle_pw = 15;
end
radius2origin = 0.6;
angle2origin  = [0 180];
%weight = 0.05;

loudspeakers   = loudspeaker_setup;  % Number of loudspeakers
if loudspeaker_setup == 16 || loudspeaker_setup == 32
    speaker_arc    = 180;  % Degrees
    first_speaker  = 90; % Degrees
else
    speaker_arc    = 360;  % Degrees
    first_speaker  = 0; % Degrees
end
speaker_radius = 1.5; % Metres

Drive = 'Z:\';
Output_file_path     = [Drive '+Speaker_Signals\']; % Can be relative or exact
Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
SetupInfo            = ['_' num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz_' ...
    num2str(angle_pw) 'pwAngle_' ...
    num2str(Noise_Mask_dB) 'dB_' ...
    num2str(weight) 'weight__with' masker_type];
Output_file_name     = [Input_file_name '__' ...
    num2str(loudspeakers) 'spkrs_' ...
    SetupInfo];
Output_file_ext      = '.WAV';





% Read input signal
Input_Signal = audioread( [Input_file_path Input_file_name Input_file_ext] );

for sig = 1:3 % Firstly we compute the loudspeaker signals for the input signal then we compute the loudspeaker signals for the additive zone weighted noise for the stereo pair
    
    %% First, Load the relevant look-up tables and check compatability
    if sig == 1
        load([Drive '+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat']);
    elseif sig == 2
        load([Drive '+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(leakage_angle_L) 'deg_' LUT_resolution ...
            '_zones_swapped.mat']);
    elseif sig == 3
        load([Drive '+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(leakage_angle_R) 'deg_' LUT_resolution ...
            '_zones_swapped.mat']);
    end
    
    %Loudspeaker capable LUT available?
    if ~exist('Loudspeaker_Weights__Weight_Vs_Frequency','var')
        error('A Look-Up Table with valid Loudspeaker Weights was not found. Please either choose another LUT or generate a valid LUT.');
    end
        
    %% Second, Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
    % The amount of aliasing is predicted from the average magnitude in the
    % quiet zone. Where there is a large amount of aliasing and hence a
    % large magnitude in the quiet zone, we invert this level and apply it
    % to the noise input signal.
    if (sig == 2) || (sig == 3)
        W = -mag2db(Quiet_Sample__Weight_Vs_Frequency);        
        freqs = linspace(0, Fs/2, length(Input_Signal)/2 + 1);
        freqs = freqs(freqs>=min(Frequencies) & freqs<=max(Frequencies));
        W_ = permute( Tools.interpVal_2D(W, Frequencies, Weights, freqs, weight, 'spline'), [2 1]);        
        
        % Equalise signal in target "bright" zone
        Input_Signal = applyWeight(Input_Signal, W_, freqs, Fs);       
        
        %Find cutoff frequencies for band pass filter
        Alias_leakage_threshold = 7.5; %dB
        [~,bandpass_centre] = max(W_);
        f_cutoff_low  = freqs(find( W_(1:bandpass_centre) <=Alias_leakage_threshold,1,'last'));
        f_cutoff_high = freqs(find( W_(bandpass_centre:end) <=Alias_leakage_threshold,1,'first'));
        f_cutoff = [f_cutoff_low, f_cutoff_high];        
        
        %Design low pass filter
        [b,a] = butter(6, f_cutoff./(Fs/2));
        
        %Apply low pass filter to noise
        Input_Signal = filter(b,a,Input_Signal(:));
    end
    
    
    %% Third, find the frequency domain representation of the audio file that is wished to be reproduced in the spatial domain.
    [Z, Frequencies_, ~, Windows] = Broadband_Tools.FFT_custom(Input_Signal, Nfft, Fs, overlap);
    
    % Truncate to frequencies in the range f_low <-> f_high
    trunc_index_low  = find(Frequencies_ < f_low , 1, 'last' ) + 1;
    trunc_index_high = find(Frequencies_ > f_high, 1 ) + 1;
    if isempty(trunc_index_low)
        trunc_index_low = 1;
    end
    if isempty(trunc_index_high)
        trunc_index_high = length(Frequencies_);
    end
    
    Frequencies_ = Frequencies_( :, trunc_index_low:trunc_index_high );
    
    
    %% Fourth, build a flat spectra desired multizone soundfield for all frequencies from the previous fft and save the speaker weights for each frequency bin.
    [szW, szF] = size(Loudspeaker_Weights__Weight_Vs_Frequency);
    LUT_Loudspeaker_Weights = cell2mat(Loudspeaker_Weights__Weight_Vs_Frequency);
    LUT_Loudspeaker_Weights = permute( reshape(LUT_Loudspeaker_Weights, loudspeakers, szW, szF), [2 3 1] );
    
    % Find the wrights that will give us the biggest contrast possible
    % (works better at lower frequencies)
    LUT_MagDiff = Bright_Sample__Weight_Vs_Frequency - Quiet_Sample__Weight_Vs_Frequency;
    Ideal_val = max(LUT_MagDiff(:))*ones(1,length(Frequencies_)); %Ideal value is the maximum separation
    weights = Tools.interpFromVal_2D(LUT_MagDiff, Frequencies, Weights, Frequencies_, Ideal_val);
            
    % When interpolating the angle of the complex loudspeaker weight we need to phase unwrap otherwise
    % the interpolation may become close to 180 degrees out of phase which will
    % cause contructive interference instead of destructive and vise versa
    Loudspeaker_Weights = zeros(length(Frequencies_),loudspeakers);
    for spkr = 1:loudspeakers
        Loudspeaker_Weights(:,spkr) = permute( Tools.interpVal_2D(LUT_Loudspeaker_Weights(:,:,spkr), Frequencies, Weights, Frequencies_, weights, 'spline'), [2 1]);
        
        Loudspeaker_Weights(:,spkr) = permute( Tools.interpVal_2D(abs(LUT_Loudspeaker_Weights(:,:,spkr)), Frequencies, Weights, Frequencies_, weights, 'spline'), [2 1]) ...
            .* exp(1i * angle(Loudspeaker_Weights(:,spkr)));
        
    end
    

    
    %% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal
    % % Here we want to build the speaker signals for each speaker so that our loudspeaker weights are taken into account.
    % % We want to form the entire spectrum by adding the conjugate of the frame
    % % to the existing frame where the negative frequencies of the transform
    % % would usually exist.
 Loudspeaker_Weights = [zeros(trunc_index_low-1, loudspeakers); ...
                        Loudspeaker_Weights; ...
                        zeros( size(Z,2) - trunc_index_high, loudspeakers)];                    
 Loudspeaker_Weights = permute( repmat(Loudspeaker_Weights, [1 1 size(Z,1)]), [3 1 2]);
 Z_l = repmat(Z, [1 1 loudspeakers]);
 
    %
    Loudspeakers_ = zeros( [size(Z_l,1) (size(Z_l,2))*2 size(Z_l,3)] );
    for spkr = 1:loudspeakers
        Loudspeakers_(:,1:end/2,spkr) = Z_l(:,:,spkr) .* Loudspeaker_Weights(:,:,spkr);
        Loudspeakers_(:,end/2+1:end,spkr) = conj( [-Loudspeakers_(:,1,spkr).*0 Loudspeakers_(:,end/2:-1:2,spkr)]);
    end
    Original = [Z(:,:,1) conj( [-Z(:,1,1).*0 Z(:,end:-1:2,1)] )];
    % clear Loudspeaker_Weights; % Save on memory
    %
    % % We then want to perform an Inverse FFT (ifft) on each full spectrum frame
    for frame = 1:size(Loudspeakers_, 1)
        for spkr = 1:loudspeakers
            Loudspeakers_(frame,:,spkr) = ifft( Loudspeakers_(frame,:,spkr) );
        end
        Original(frame,:) = ifft( Original(frame, :) );
    end
    
    %We should apply the second square root hamming window here
    %we should do this to remove artificats caused by our spectral
    %modification
    %for frame = 1:size(Loudspeakers_, 1)
    for spkr = 1:loudspeakers
        Loudspeakers_(:,:,spkr) = Loudspeakers_(:,:,spkr) .* Windows;
    end
    Original = Original .* Windows;
    %end
    
    %
    % % Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
    % %Loudspeaker_Signals =
    % zeros([(size(Z,1)+ceil(overlap))*size(Z,2)*2*(1-overlap) loudspeakers] ); % pre-allocate memory
    for spkr = 1:loudspeakers
        Loudspeaker_Signals(:,spkr) = Broadband_Tools.OverlapAdd( Loudspeakers_(:,:,spkr), overlap ); %#ok<AGROW>
    end
    Original_ = Broadband_Tools.OverlapAdd( Original, overlap );
    % clear Loudspeakers_; % Save on memory
    
    
    % Scale signals so they don't clip upon saving
    if Noise_Mask_dB <= 0
        scaler = 1 / (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
    elseif Noise_Mask_dB > 0 % For a positive masker we scale the signals to save up to a 40db noise masker
        scaler = 1 / (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
    end
    
    % Normalise Loudspeaker Signals
    Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)))  *  scaler; 
    Original_ = Original_ ./ max(abs(Original_(:)))  *  scaler;
    
    if sig == 1
        Loudspeaker_Signals_Input_Signal = Loudspeaker_Signals;
        Original_Input_Signal = Original_;
    elseif sig == 2
        Loudspeaker_Signals_Noise_L = Loudspeaker_Signals;
    elseif sig == 3
        Loudspeaker_Signals_Noise_R = Loudspeaker_Signals;
    end
    
    
    max_Spkrval =  max( abs( Loudspeaker_Signals_Input_Signal(:) ) );
    level_mag = db2mag(Noise_Mask_dB);
    Input_Signal = Perceptual_Tools.GreyNoise( length(Loudspeaker_Signals_Input_Signal)/Fs, Fs, max_Spkrval * level_mag );
    
end
Original_ = Original_Input_Signal;


%% Add Stereo Zone Weighted Noise Loudspeaker Signals to Speech Loudspeaker Signals
 % Here we add our loudspeaker signals which reproduce our input signal and a relatively weighted stereo 'zone weighted' noise
 % That is to say, we add the stereo zone weighted noise to the reproduction,
 % however, the stereo zone weighted noise can be varied in level such that 0dB is
 % when the peak noise is equal to the input signal reproduction.
Loudspeaker_Signals = Loudspeaker_Signals_Input_Signal + db2mag(Noise_Mask_dB) * (Loudspeaker_Signals_Noise_L + Loudspeaker_Signals_Noise_R) / 2;



%% Once we have the speaker signals we should save them for later use as .wav files
filenumbers = num2str((1:loudspeakers)');
filenumbers(filenumbers==' ') = '_';
fullpath = [repmat([Output_file_path Output_file_path_ext Output_file_name], [loudspeakers 1]) ...
    filenumbers ...
    repmat(Output_file_ext, [loudspeakers 1]) ];

if ~exist([Output_file_path Output_file_path_ext],'dir'); mkdir([Output_file_path Output_file_path_ext]); end

for spkr = 1:loudspeakers
    audiowrite(fullpath(spkr,:), Loudspeaker_Signals(:, spkr), Fs);
end
audiowrite([Output_file_path Output_file_path_ext ...
    Input_file_name '__Original_' ...
    num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz' ...
    Input_file_ext], Original_, Fs);




end

function y = applyWeight(x, W, W_freqs, Fs)

% Frequency Domain Tranform
X = fft(x); % Fast Fourier Transform

% Frequencies
M = length(x);
NumPts = M/2 + 1;
freqs = linspace(0, Fs/2, NumPts);
cutoff_low = min(W_freqs);
cutoff_high = max(W_freqs);

% Weight Levels
W = [linspace(0, W(1), length(freqs(freqs<cutoff_low))), ...
     W', ...
     linspace(W(end), 0, length(freqs(freqs>cutoff_high)))];

% Apply magnitude weighting
X(1:NumPts) = X(1:NumPts) .* db2mag(W);

% Apply conjugation for negative frequency side of spectrum
X(NumPts+1:M) = conj(X(M/2:-1:2));

% Time Domain Transform
y = ifft(X); % Inverse Fast Fourier Transform

% prepare output vector y
y = real(y(1, 1:M));

% remove DC
y = y(:) - mean(y);
end