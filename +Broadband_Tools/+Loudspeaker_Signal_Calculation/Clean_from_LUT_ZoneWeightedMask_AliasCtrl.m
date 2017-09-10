function Clean_from_LUT_ZoneWeightedMask_AliasCtrl( Input_file, LUT_resolution, Noise_Mask_dB, weight, setup, setup_mask )
%Clean_from_LUT_ZoneWeightedMask Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables
Drive = 'Z:\';
[Input_file_path, Input_file_name, Input_file_ext] = fileparts( Input_file );
Input_file_path = [Input_file_path '\'];

signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.L_noise_mask = Noise_Mask_dB; % dB
signal_info.weight = weight;
signal_info.method = 'ZoneWeightMaskAliasCtrl';
signal_info.input_filename = Input_file_name;

[Output_path, Output_file_name, Output_file_ext] = ...
    Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, LUT_resolution, Drive, 'new');


loudspeakers   = setup.Loudspeaker_Count;



% Read input signal
Input_Signal = audioread( Input_file );

for sig = 1:2 % Firstly we compute the loudspeaker signals for the input signal then we compute the loudspeaker signals for the additive zone weighted noise
    
    %% First, Load the relevant look-up tables and check compatability
    method = {'new4', 'new3', 'new2', 'new'};
    for m = 1:2
        if sig == 1
            method_ = [method, {'old'}];
            [DB,err] = Soundfield_Database.loadDatabaseFromSetup( setup, LUT_resolution, Drive, method_{m} );
        elseif sig == 2
            method_ = [method, {'old_zones_swapped'}];
            [DB,err] = Soundfield_Database.loadDatabaseFromSetup( setup_mask, LUT_resolution, Drive, method_{m} );
        end
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
    noise_freqs = linspace(0, Fs/2, len/2 + 1);
    noise_freqs = noise_freqs(noise_freqs>=min(Frequencies) & noise_freqs<=max(Frequencies));
    len = Nfft;
    freqs = linspace(0, Fs/2, len/2 + 1);
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
    if sig == 2
%         W = -mag2db(DB.Quiet_Sample__Weight_Vs_Frequency);   
%         W_ = permute( Tools.interpVal_2D(W, Frequencies, Weights, noise_freqs, noise_weights, 'spline'), [2 1]);        
%         W_=W_(:);
        
        % Equalise signal in target "bright" zone
        %Input_Signal = applyWeight(Input_Signal, W_, noise_freqs, Fs);       
        
        %Find cutoff frequencies for band pass filter
%         Alias_leakage_threshold = 7.5; %dB
%         [~,bandpass_centre] = max(W_);
%         f_cutoff_low  = noise_freqs(find( W_(1:bandpass_centre) <=Alias_leakage_threshold,1,'last'));
%         f_cutoff_high = noise_freqs( (bandpass_centre-1) + find( W_(bandpass_centre:end) <=Alias_leakage_threshold,1,'first'));
%         f_cutoff = [f_cutoff_low, f_cutoff_high];
        R_ = max( [setup.Multizone_Soundfield.Bright_Zone.Radius_q + setup.Multizone_Soundfield.Bright_Zone.Origin_q.Distance; ...
                   setup.Multizone_Soundfield.Quiet_Zone.Radius_q + setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Distance;  ]);
        phiL_rad = setup.Speaker_Arc_Angle / 180 * pi;
        
        f_cutoff = c * (loudspeakers - 1) / (2 * R_ * phiL_rad);
        
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
        Loudspeaker_Signals(:,spkr) = Tools.OverlapAdd( Loudspeakers_(:,:,spkr), overlap ); %#ok<AGROW>
    end
    Original_ = Tools.OverlapAdd( Original, overlap );
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
    end
    
    %Input_Signal = Tools.generateNoise(Loudspeaker_Signals, Noise_Mask_dB, 'WGN'); % Generate noise to put back into the system and zone weight according to the multizone setup
    max_Spkrval =  max( abs( Loudspeaker_Signals(:) ) );
    level_mag = db2mag(Noise_Mask_dB);
    Input_Signal = Perceptual_Tools.GreyNoise( length(Loudspeaker_Signals)/Fs, Fs, max_Spkrval * level_mag );
    Input_Signal = Input_Signal(:);
end
Original_ = Original_Input_Signal;


%% Add Zone Weighted Noise Loudspeaker Signals to Speech Loudspeaker Signals
 % Here we add our loudspeaker signals which reproduce our input signal and a relatively weighted 'zone weighted' noise
 % That is to say, we add the zone weighted noise to the reproduction,
 % however, the zone weighted noise can be varied in level such that 0dB is
 % when the peak noise is equal to the input signal reproduction.
Loudspeaker_Signals = Loudspeaker_Signals_Input_Signal + db2mag(Noise_Mask_dB) * Loudspeaker_Signals;






%% Once we have the speaker signals we should save them for later use as .wav files
filenumbers = num2str((1:loudspeakers)');
filenumbers(filenumbers==' ') = '_';
fullpath = [repmat([Output_path Output_file_name], [loudspeakers 1]) ...
    filenumbers ...
    repmat(Output_file_ext, [loudspeakers 1]) ];

if ~exist(Output_path,'dir'); mkdir(Output_path); end

for spkr = 1:loudspeakers
    audiowrite(fullpath(spkr,:), Loudspeaker_Signals(:, spkr), Fs);
end
audiowrite([Output_path ...
    Input_file_name '_' 'Original'
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