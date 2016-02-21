function ZoneWeightedMasker_AliasCtrl( Signal_Length, LUT_resolution, Noise_Mask_dB, weight, setup )
%Clean_from_LUT_ZoneWeightedMask Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables
Drive = 'Z:\';

signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.L_noise_mask = Noise_Mask_dB; % dB
signal_info.weight = weight;
signal_info.method = 'ZoneWeightMaskerAliasCtrl';
signal_info.input_filename = 'Masker';

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, LUT_resolution, Drive, 'new');


loudspeakers   = setup.Loudspeaker_Count;



% Generate Input Signal
level_mag = db2mag(signal_info.L_noise_mask);
Input_Signal = Perceptual_Tools.GreyNoise( Signal_Length/signal_info.Fs, signal_info.Fs, level_mag );

% Compute the loudspeaker signals for the additive zone weighted noise    
    %% First, Load the relevant look-up tables and check compatability
    method = {'new4', 'new3', 'new2', 'new'};
    for m = 1:2
            method_ = [method, {'old_zones_swapped'}];
            [DB,err] = Soundfield_Database.loadDatabaseFromSetup( setup, LUT_resolution, Drive, method_{m} );
        if ~err
            break;
        end
    end
    
    %Loudspeaker capable LUT available?
    if ~isfield(DB,'Loudspeaker_Weights__Weight_Vs_Frequency')
        error('A Look-Up Table with valid Loudspeaker Weights was not found. Please either choose another LUT or generate a valid LUT.');
    end
    
    Frequencies = DB.Frequencies;
    Weights = DB.Weights;
    
    %% Find ideal weights
    single_weight = false;
    
    len = length(Input_Signal);
    noise_freqs = linspace(0, signal_info.Fs/2, len/2 + 1);
    noise_freqs = noise_freqs(noise_freqs>=min(Frequencies) & noise_freqs<=max(Frequencies));
    len = signal_info.Nfft;
    freqs = linspace(0, signal_info.Fs/2, len/2 + 1);
    freqs = freqs(freqs>=min(Frequencies) & freqs<=max(Frequencies));
    
    if single_weight
        noise_weights = repmat(weight,1,length(noise_freqs));
        weights = repmat(weight,1,length(freqs));
    else
        % Find the weights that will give us the biggest contrast possible
        % (works better at lower frequencies)
        LUT_MagDiff = DB.Acoustic_Contrast__Weight_Vs_Frequency;%DB.Bright_Sample__Weight_Vs_Frequency - DB.Quiet_Sample__Weight_Vs_Frequency;
        
        % Noise weights
        LUT_MagDiff_interp = interp2(Frequencies,Weights,LUT_MagDiff,noise_freqs',Weights,'spline');
        [~,I]=max(LUT_MagDiff_interp);
        noise_weights = Weights(I);
        
        %Signal Weights
        LUT_MagDiff_interp = interp2(Frequencies,Weights,LUT_MagDiff,freqs',Weights,'spline');
        [~,I]=max(LUT_MagDiff_interp);
        weights = Weights(I);
    end
    
    %% Second, Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
    % The amount of aliasing is predicted from the average magnitude in the
    % quiet zone. Where there is a large amount of aliasing and hence a
    % large magnitude in the quiet zone, we invert this level and apply it
    % to the noise input signal.
    
%         W = -mag2db(DB.Quiet_Sample__Weight_Vs_Frequency);   
%         W_ = permute( Tools.interpVal_2D(W, Frequencies, Weights, noise_freqs, noise_weights, 'spline'), [2 1]);        
%         W_=W_(:);
        
        % Equalise signal in target "bright" zone
        %Input_Signal = applyWeight(Input_Signal, W_, noise_freqs, signal_info.Fs);       
        
        %Find cutoff frequencies for band pass filter
%         Alias_leakage_threshold = 7.5; %dB
%         [~,bandpass_centre] = max(W_);
%         f_cutoff_low  = noise_freqs(find( W_(1:bandpass_centre) <=Alias_leakage_threshold,1,'last'));
%         f_cutoff_high = noise_freqs( (bandpass_centre-1) + find( W_(bandpass_centre:end) <=Alias_leakage_threshold,1,'first'));
%         f_cutoff = [f_cutoff_low, f_cutoff_high];

        % This method formulates the cutoff frequency where aliasing begins to
        % occur. It is based on the frequency where aliasing occurs for the
        % reproduction region but not necessarily the individual zones.

        R_ = max( [setup.Multizone_Soundfield.Bright_Zone.Radius_q + setup.Multizone_Soundfield.Bright_Zone.Origin_q.Distance; ...
                   setup.Multizone_Soundfield.Quiet_Zone.Radius_q + setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Distance;  ]);
        phiL_rad = setup.Speaker_Arc_Angle / 180 * pi;
        
        f_cutoff = signal_info.c * (loudspeakers - 1) / (2 * R_ * phiL_rad);
        
        %Design LOW pass filter
        [b,a] = butter(6, f_cutoff./(signal_info.Fs/2), 'low');
        
        %Apply low pass filter to noise
        Input_Signal = filter(b,a,Input_Signal(:));

    
    
    %% Third, find the frequency domain representation of the audio file that is wished to be reproduced in the spatial domain.
    [Z, Frequencies_, ~, Windows] = Broadband_Tools.FFT_custom(Input_Signal, signal_info.Nfft, signal_info.Fs, signal_info.overlap);
    
    % Truncate to frequencies in the range f_low <-> f_high
    trunc_index_low  = find(Frequencies_ < signal_info.f_low , 1, 'last' ) + 1;
    trunc_index_high = find(Frequencies_ > signal_info.f_high, 1 ) + 1;
    if isempty(trunc_index_low)
        trunc_index_low = 1;
    end
    if isempty(trunc_index_high)
        trunc_index_high = length(Frequencies_);
    end
    
    Frequencies_ = Frequencies_( :, trunc_index_low:trunc_index_high );
    
%     Frequencies_ = Frequencies_(Frequencies_>=min(Frequencies) & Frequencies_<=max(Frequencies));
    
    %% Fourth, build a flat spectra desired multizone soundfield for all frequencies from the previous fft and save the speaker weights for each frequency bin.
    [szW, szF] = size(DB.Loudspeaker_Weights__Weight_Vs_Frequency);
    LUT_Loudspeaker_Weights = cell2mat(DB.Loudspeaker_Weights__Weight_Vs_Frequency);
    LUT_Loudspeaker_Weights = permute( reshape(LUT_Loudspeaker_Weights, loudspeakers, szW, szF), [2 3 1] );    
    
    % When interpolating the angle of the complex loudspeaker weight we need to phase unwrap otherwise
    % the interpolation may become close to 180 degrees out of phase which will
    % cause contructive interference instead of destructive and vise versa
    Loudspeaker_Weights = zeros(length(Frequencies_),loudspeakers);
    for spkr = 1:loudspeakers
        LW = Tools.interpVal_2D(LUT_Loudspeaker_Weights(:,:,spkr), Frequencies, Weights, Frequencies_, weights, 'spline');
        Loudspeaker_Weights(:,spkr) = LW(:);
        
        LW_abs = Tools.interpVal_2D( abs(LUT_Loudspeaker_Weights(:,:,spkr)), Frequencies, Weights, Frequencies_, weights, 'spline');
        Loudspeaker_Weights(:,spkr) = LW_abs(:) ...
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
        Loudspeaker_Signals(:,spkr) = Broadband_Tools.OverlapAdd( Loudspeakers_(:,:,spkr), signal_info.overlap ); %#ok<AGROW>
    end
    Original_ = Broadband_Tools.OverlapAdd( Original, signal_info.overlap );
    % clear Loudspeakers_; % Save on memory
    
    
    % Scale signals so they don't clip upon saving
    if signal_info.L_noise_mask <= 0
        scaler = 1 / (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
    elseif signal_info.L_noise_mask > 0 % For a positive masker we scale the signals to save up to a 40db noise masker
        scaler = 1 / (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
    end
    
    % Normalise Loudspeaker Signals
    Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)))  *  scaler  *  level_mag; 
    Original_ = Original_ ./ max(abs(Original_(:)))  *  scaler;
    





%% Once we have the speaker signals we should save them for later use as .wav files
Broadband_Tools.Loudspeaker_Signal_Calculation.saveLoudspeakerSignals( ...
    Output_path, ...
    Output_file_name, ...
    Output_file_ext, ...
    Loudspeaker_Signals, ...
    [], [], [], [], ...
    signal_info.Fs );


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
     W(:)', ...
     linspace(W(end), 0, length(freqs(freqs>cutoff_high)))];

% Apply magnitude weighting
X(1:NumPts) = X(1:NumPts) .* db2mag(W(:));

% Apply conjugation for negative frequency side of spectrum
X(NumPts+1:M) = conj(X(M/2:-1:2));

% Time Domain Transform
y = ifft(X); % Inverse Fast Fourier Transform

% prepare output vector y
y = real(y(1:M));

% remove DC
y = y(:) - mean(y);
end