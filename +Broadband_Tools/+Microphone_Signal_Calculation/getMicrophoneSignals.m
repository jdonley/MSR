function [ Mic_Signals ] = getMicrophoneSignals( Input_Signal, SYS )
% Summary of this function goes here
% 
% Syntax:	[ Mic_Signals, Original_ ] = getMicrophoneSignals( Input_Signal, SYS )
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% Outputs: 
% 	output1 - Description
% 	output2 - Description
% 
% Example: 
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 17 August 2017
% Version: 0.1 (17 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB = Room_Acoustics.loadRIRDatabaseFromSetup(SYS);

Mic_Signals = Tools.fconv( ...
    repmat(Input_Signal,1,SYS.Room_Setup.NoReceivers), ...
    DB.RIRs.Bright_RIRs.');

%{
if ~isfield(signal_info,'f_high_meas')
    signal_info.f_high_meas = signal_info.f_high;
end
if ~isfield(signal_info,'f_low_meas')
    signal_info.f_low_meas = signal_info.f_low;
end

%% First, Load the relevant look-up tables and check compatability
method = {'new4', 'new3', 'new2', 'new'};
for m = 1:2
    method_ = [method, {'old_zones_swapped'}];
    [DB,err] = Soundfield_Database.loadDatabaseFromSetup( setup, system_info.LUT_resolution, system_info.Drive, method_{m} );
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

len = length(Input_Signal);
noise_freqs = linspace(0, signal_info.Fs/2, len/2 + 1);
noise_freqs = noise_freqs(noise_freqs>=min(Frequencies) & noise_freqs<=max(Frequencies));

len = signal_info.Nfft + signal_info.zeropadtime * signal_info.Fs;
freqs = linspace(0, signal_info.Fs/2, len/2 + 1);
freqs = freqs(freqs>=min(Frequencies) & freqs<=max(Frequencies));

if ~strcmpi(signal_info.weight, 'auto')
    noise_weights = repmat(signal_info.weight,1,length(noise_freqs));
    weights = repmat(signal_info.weight,1,length(freqs));
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


%% Second, find the frequency domain representation of the audio file that is wished to be reproduced in the spatial domain.
si = signal_info;
if isempty(strfind(lower(signal_info.method),'cancel'))
    si.time_delay = 0;
end
[Z, Frqs, ~, Windows, Zo] = Broadband_Tools.FFT_custom( Input_Signal, si);

% Truncate to frequencies in the range f_low <-> f_high
iLo  = find(Frqs < signal_info.f_low , 1, 'last' ) + 1;
iHi = find(Frqs > signal_info.f_high, 1 ) + 1;
if isempty(iLo)
    iLo = 1;
end
if isempty(iHi)
    iHi = length(Frqs);
end

FI = iLo:iHi;
Frequencies_ = Frqs( :, FI );

%     Frequencies_ = Frequencies_(Frequencies_>=min(Frequencies) & Frequencies_<=max(Frequencies));

%% Third, build a flat spectra desired multizone soundfield for all frequencies from the previous fft and save the speaker weights for each frequency bin.
[szW, szF] = size(DB.Loudspeaker_Weights__Weight_Vs_Frequency);
LUT_Loudspeaker_Weights = cell2mat(DB.Loudspeaker_Weights__Weight_Vs_Frequency);
LUT_Loudspeaker_Weights = permute( reshape(LUT_Loudspeaker_Weights, setup.Loudspeaker_Count, szW, szF), [2 3 1] );

% When interpolating the angle of the complex loudspeaker weight we need to phase unwrap otherwise
% the interpolation may become close to 180 degrees out of phase which will
% cause contructive interference instead of destructive and vise versa
Loudspeaker_Weights = zeros(length(Frequencies_),setup.Loudspeaker_Count);
ind = 1:(iHi-iLo+1);
for spkr = 1:setup.Loudspeaker_Count
    LW = Tools.interpVal_2D(LUT_Loudspeaker_Weights(:,ind,spkr), Frequencies(ind), Weights, Frequencies_, weights, 'spline');
    Loudspeaker_Weights(:,spkr) = LW(:);
    
    LW_abs = Tools.interpVal_2D( abs(LUT_Loudspeaker_Weights(:,ind,spkr)), Frequencies(ind), Weights, Frequencies_, weights, 'spline');
    Loudspeaker_Weights(:,spkr) = LW_abs(:) ...
        .* exp(1i * angle(Loudspeaker_Weights(:,spkr)));
end


%% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal
% % Here we want to build the speaker signals for each speaker so that our loudspeaker weights are taken into account.
% % We want to form the entire spectrum by adding the conjugate of the frame
% % to the existing frame where the negative frequencies of the transform
% % would usually exist.
N = signal_info.Nfft;

taperWin = tukeywin(N,(iLo+0.05*N)/N);

% test_Weights = [Loudspeaker_Weights(1,8).*ones(iLo-1, 1); ... %For bandpass filtering the original to maintain fair comparison
%     Loudspeaker_Weights(:,8); ...
%     zeros( N/2 - iHi, 1)];

Loudspeaker_Weights = [repmat(Loudspeaker_Weights(1,:),iLo-1,1); ...
    Loudspeaker_Weights;...
    zeros( N/2 - iHi, setup.Loudspeaker_Count)];
Orig_Weights = [ones(iLo-1, 1); ... %For bandpass filtering the original to maintain fair comparison
    ones(length(Frequencies_),1); ...
    zeros( N/2 - iHi, 1)];

% test_Weights = test_Weights .* taperWin;
% Loudspeaker_Weights = Loudspeaker_Weights .* repmat(taperWin,1,setup.Loudspeaker_Count);
% Orig_Weights = Orig_Weights .* taperWin;

Loudspeaker_Weights = permute( repmat(Loudspeaker_Weights, [1 1 size(Z,1)]), [3 1 2]);
Orig_Weights = permute( repmat(Orig_Weights, [1 1 size(Zo,1)]), [3 1 2]);
% test_Weights = permute( repmat(test_Weights, [1 1 size(Z,1)]), [3 1 2]);
Z_l = repmat(Z, [1 1 setup.Loudspeaker_Count]);

%
Loudspeakers_ = zeros( [size(Z_l,1) (size(Z_l,2))*2 size(Z_l,3)] );
for spkr = 1:setup.Loudspeaker_Count
    Loudspeakers_(:,1:end/2,spkr) = Z_l(:,:,spkr) .* Loudspeaker_Weights(:,:,spkr);
    Loudspeakers_(:,end/2+1:end,spkr) = conj( [-Loudspeakers_(:,1,spkr).*0 Loudspeakers_(:,end/2:-1:2,spkr)]);
end
Original = [Zo(:,:,1) conj( [-Zo(:,1,1).*0 Zo(:,end:-1:2,1)] )] .* [Orig_Weights Orig_Weights(:,1) Orig_Weights(:,end:-1:2)];
%  Input_toMatch = [Z(:,:,1) conj( [-Z(:,1,1).*0 Z(:,end:-1:2,1)] )].* [test_Weights test_Weights(:,1) test_Weights(:,end:-1:2)];

%
% % We then want to perform an Inverse FFT (ifft) on each full spectrum frame
for frame = 1:size(Loudspeakers_, 1)
    for spkr = 1:setup.Loudspeaker_Count
        Loudspeakers_(frame,:,spkr) = ifft( Loudspeakers_(frame,:,spkr) * size(Loudspeakers_,2) );
    end
end
for frame = 1:size(Original,1)
    Original(frame,:) = ifft( Original(frame, :) * size(Original,2) );
end
% for frame = 1:size(Input_toMatch,1)
%     Input_toMatch(frame,:) = ifft( Input_toMatch(frame, :) * size(Input_toMatch,2) );
% end


%We should apply the second square root hamming window here
%we should do this to remove artificats caused by our spectral
%modification
%for frame = 1:size(Loudspeakers_, 1)
for spkr = 1:setup.Loudspeaker_Count
    Loudspeakers_(:,:,spkr) = Loudspeakers_(:,:,spkr) .* Windows;
end
Original = Original .* Windows(1:size(Original,1),:);
% Input_toMatch = Input_toMatch .* Windows(1:size(Input_toMatch,1),:);
%end

if ~isempty(signal_info.zeropadtime)
    Npad=signal_info.zeropadtime * signal_info.Fs;
    Loudspeakers_(:,[1:Npad/2, end-Npad/2+1:end],:)=[];
    Original(:,[1:Npad/2, end-Npad/2+1:end],:)=[];
    %     Input_toMatch(:,[1:Npad/2, end-Npad/2+1:end],:)=[];
end

%
% % Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
% %Loudspeaker_Signals =
% zeros([(size(Z,1)+ceil(overlap))*size(Z,2)*2*(1-overlap) setup.Loudspeaker_Count] ); % pre-allocate memory
for spkr = 1:setup.Loudspeaker_Count
    %     Loudspeaker_Signals(:,spkr) = Broadband_Tools.OverlapAdd( real(Loudspeakers_(:,:,spkr)), signal_info.overlap ); %#ok<AGROW>
    Mic_Signals(:,spkr) = overlapadd( squeeze(real(Loudspeakers_(:,:,spkr))), ones(signal_info.Nfft,1), (1-signal_info.overlap)*signal_info.Nfft  ); %#ok<AGROW>
    Mic_Signals(isnan(Mic_Signals(:,spkr)),spkr)=0;
end
% Original_ = Broadband_Tools.OverlapAdd( Original, signal_info.overlap );
Original_ = overlapadd( Original, ones(signal_info.Nfft,1), (1-signal_info.overlap)*signal_info.Nfft  ); %#ok<AGROW>
% Input_toMatch_ = overlapadd( Input_toMatch, ones(signal_info.Nfft,1), (1-signal_info.overlap)*signal_info.Nfft  ); %#ok<AGROW>
% Input_toMatch_ = Broadband_Tools.OverlapAdd( Input_toMatch, signal_info.overlap );
% clear Loudspeakers_; % Save on memory



%%
% if ~isempty(signal_info.time_delay)
%     t_d = signal_info.time_delay;
% else
%     t_d = signal_info.Nfft/signal_info.Fs;
% end
%     t_d = signal_info.time_delay;
% if t_d ~= 0 && ~isempty(strfind(lower(signal_info.method),'cancel'))
% if ~isempty(strfind(lower(signal_info.method),'cancel'))
if signal_info.predict_buff ~= 0
    LS_=[];
    for s=1:setup.Loudspeaker_Count
        LS = enframe(Mic_Signals(:,s),signal_info.Nfft/2,signal_info.Nfft/2);
        LS_(:,s) = reshape(LS(2:2:end,:).',[],1);
    end
    Mic_Signals = LS_;
    Mic_Signals(size(Mic_Signals,1):length(Original_),:)=0;
end
%     IM = enframe(Input_toMatch_,signal_info.Nfft/2,signal_info.Nfft/2);
%     IM_ = reshape(IM(2:2:end,:).',[],1);
%     Input_toMatch_ = IM_;
%     Input_toMatch_(size(Input_toMatch_,1):length(Original_),:)=0;
% end
% s2=size(Loudspeaker_Signals,2);
% Zs = zeros(t_d*signal_info.Fs,1);
% if ~isempty(strfind(lower(signal_info.method),'cancel'))
%     Loudspeaker_Signals = [repmat(Zs,1,s2); ...
%         Loudspeaker_Signals];
%     Original_ = [Original_; Zs];
% else
%     Loudspeaker_Signals = [Loudspeaker_Signals; ...
%         repmat(Zs,1,s2)];
%     Original_ = [Original_; Zs];
% end

%%

% Scale signals so they don't clip upon saving
if signal_info.L_noise_mask <= 0
    scaler = 1 / (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
elseif signal_info.L_noise_mask > 0 % For a positive masker we scale the signals to save up to a 40db noise masker
    scaler = 1 / (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
end

% Normalise Loudspeaker Signals
% [~,adjVal] = Broadband_Tools.power_norm( Input_Signal(:), Input_toMatch_(:), signal_info.Fs, [signal_info.f_low_meas, signal_info.f_high_meas]);
Mic_Signals = Mic_Signals  *  scaler;

%}
end

