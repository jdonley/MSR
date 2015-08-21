function Clean_from_LUT_ZoneWeightedMask( Input_file_path, Input_file_name, Input_file_ext, LUT_resolution, Noise_Mask_dB, weight, loudspeaker_setup, angle_pw )
%Clean_from_LUT_ZoneWeightedMask Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables

Fs = 16000; % Sampling frequency
Nfft = 1024;% Number of fft components
overlap = 0.5;
f_low  = 150;  % Hz
f_high = 8000; % Hz

resolution = 100;
phase = 0;
radius = 0.3;
if nargin < 8
    angle_pw = 15;
end
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

Drive = 'Z:\';
Output_file_path     = [Drive '+Speaker_Signals\']; % Can be relative or exact
Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
SetupInfo            = ['_' num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz_' ...
    num2str(angle_pw) 'pwAngle_' ...
    num2str(Noise_Mask_dB) 'dB_' ...
    num2str(weight) 'weight__withFlatMask'];
Output_file_name     = [Input_file_name '__' ...
    num2str(loudspeakers) 'spkrs_' ...
    SetupInfo];
Output_file_ext      = '.WAV';










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

%Loudspeaker capable LUT available?
if ~exist('Loudspeaker_Weights__Weight_Vs_Frequency','var')
    error('A Look-Up Table with valid Loudspeaker Weights was not found. Please either choose another LUT or generate a valid LUT.');
end

LUT_Bright = Bright_Sample__Weight_Vs_Frequency;%(db2mag(Contrast__Weight_Vs_Frequency).^0.5 .* db2mag(Quiet_SPL__Weight_Vs_Frequency-94));
LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;%db2mag(Quiet_SPL__Weight_Vs_Frequency-94);

[szW, szF] = size(Loudspeaker_Weights__Weight_Vs_Frequency);
LUT_Loudspeaker_Weights = cell2mat(Loudspeaker_Weights__Weight_Vs_Frequency);
LUT_Loudspeaker_Weights = permute( reshape(LUT_Loudspeaker_Weights, loudspeakers, szW, szF), [2 3 1] );

Bright_Zone_sample = permute( Tools.interpVal_2D(LUT_Bright, Frequencies, Weights, Frequencies_, weight), [2 1]);
Quiet_Zone_sample  = permute( Tools.interpVal_2D(LUT_Quiet, Frequencies, Weights, Frequencies_, weight), [2 1]);

% When interpolating the angle of the complex loudspeaker weight we need to phase unwrap otherwise
% the interpolation may become close to 180 degrees out of phase which will
% cause contructive interference instead of destructive and vise versa
Loudspeaker_Weights = zeros(length(Frequencies_),loudspeakers);
for spkr = 1:loudspeakers
    Loudspeaker_Weights(:,spkr) = permute( Tools.interpVal_2D(LUT_Loudspeaker_Weights(:,:,spkr), Frequencies, Weights, Frequencies_, weight, 'spline'), [2 1]);
    
       
%     figure(1);
%     plot(Frequencies, abs(LUT_Loudspeaker_Weights(end,:,spkr))); hold on;
%     plot(Frequencies_, abs(Loudspeaker_Weights(:,spkr))); hold off;
%     figure(2);
%     plot(Frequencies, unwrap(angle(LUT_Loudspeaker_Weights(end,:,spkr)))); hold on;
%     plot(Frequencies_, unwrap(angle(Loudspeaker_Weights(:,spkr)))); hold off;
        
    Loudspeaker_Weights(:,spkr) = permute( Tools.interpVal_2D(abs(LUT_Loudspeaker_Weights(:,:,spkr)), Frequencies, Weights, Frequencies_, weight, 'spline'), [2 1]) ...
                                .* exp(1i * angle(Loudspeaker_Weights(:,spkr)));
               
% 	figure(3);
%     plot(Frequencies, abs(LUT_Loudspeaker_Weights(end,:,spkr))); hold on;
%     plot(Frequencies_, abs(Loudspeaker_Weights(:,spkr))); hold off;
%     figure(4);
%     plot(Frequencies, unwrap(angle(LUT_Loudspeaker_Weights(end,:,spkr)))); hold on;
%     plot(Frequencies_, unwrap(angle(Loudspeaker_Weights(:,spkr)))); hold off;
                            
end


%if ~exist(['+Soundfield_Database\+From_LUT\' Output_file_path_ext],'dir'); mkdir(['+Soundfield_Database\+From_LUT\' Output_file_path_ext]); end
%save(['+Soundfield_Database\+From_LUT\' Output_file_path_ext 'Weights_and_Samples__' SetupInfo '.mat'], ...
%    'Loudspeaker_Weights', ...
%    'Bright_Zone_sample', ...
%    'Quiet_Zone_sample' );










%% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal
%load(['+Soundfield_Database\+From_LUT\' Output_file_path_ext 'Weights_and_Samples__' SetupInfo '.mat']);

% % Here we want to build the speaker signals for each speaker so that our loudspeaker weights are taken into account.
% % We want to form the entire spectrum by adding the conjugate of the frame
% % to the existing frame where the negative frequencies of the transform
% % would usually exist.
 Loudspeaker_Weights = [zeros(trunc_index_low-1, loudspeakers); Loudspeaker_Weights; zeros( size(Z,2) - trunc_index_high, loudspeakers) ];
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
% 
% % Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
% %Loudspeaker_Signals =
% zeros([(size(Z,1)+ceil(overlap))*size(Z,2)*2*(1-overlap) loudspeakers] ); % pre-allocate memory
for spkr = 1:loudspeakers
    Loudspeaker_Signals(:,spkr) = +Broadband_Tools.OverlapAdd( Loudspeakers_(:,:,spkr), overlap );
end
 Original_ = +Broadband_Tools.OverlapAdd( Original, overlap );
% clear Loudspeakers_; % Save on memory



 % Normalise Loudspeaker Signals
 Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)))  *  0.5; % The half is so that when we add 0dB of relative noise we don't clip upon saving the file
 Original_ = Original_ ./ max(abs(Original_(:)))  *  0.5;
 
 % Add white noise to Loudspeaker Signals
 Loudspeaker_Signals = Tools.addNoise(Loudspeaker_Signals, Noise_Mask_dB, 'UWGN'); 


 % % Once we have the speaker signals we should save them for later use as .wav files
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