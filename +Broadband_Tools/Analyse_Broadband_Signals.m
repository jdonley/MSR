function Analyse_Broadband_Signals( Input_file_path, Input_file_name, Input_file_ext, weight )
%FFT_CUSTOM Summary of this function goes here
%   Detailed explanation goes here



%% Initialise
%clc;
%clear;
close all;
tic;
if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))


%% Setup Variables
%Input_file_path = '+Miscellaneous\'; % Can be relative or exact
%Input_file_name = 'MaleSpeech16k';
%Input_file_ext  = '.wav';

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

Reproduction_Radius = 1.0;
f_max = 8000;
loudspeakers   = 2*ceil(f_max/343 * exp(1) * Reproduction_Radius / 2) + 1;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
first_speaker  = 0; % Degrees
speaker_radius = 1.5; % Metres

Output_file_path     = '+Speaker_Signals\'; % Can be relative or exact
Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\'];
SetupInfo            = [num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz_' ...
    num2str(angle_pw) 'pwAngle_' ...
    num2str(weight) 'weight__'];
Output_file_name     = [Input_file_name '__' ...
    num2str(loudspeakers) 'spkrs_' ...
    SetupInfo];
Output_file_ext      = '.wav';



%% Firstly, find the frequency domain representation of an audio file that is wished to be reproduced in the spatial domain.
[Z, Frequencies] = Broadband_Tools.FFT_custom([Input_file_path Input_file_name Input_file_ext], Nfft, Fs, overlap);

% Truncate to frequencies in the range f_low <-> f_high
trunc_index_low  = find(Frequencies < f_low , 1, 'last' ) + 1;
trunc_index_high = find(Frequencies > f_high, 1 ) + 1;
if isempty(trunc_index_low)
    trunc_index_low = 1;
end
if isempty(trunc_index_high)
    trunc_index_high = length(Frequencies);
end
%Z = Z ( :, trunc_index_low:trunc_index_high );
Frequencies = Frequencies( :, trunc_index_low:trunc_index_high );

% Work out the coresponding number of planewaves for these frequencies
f_samples = length(Frequencies);
Nplanewaves = round( Tools.interpVal(ceil((300-30)/((f_samples-1)^4)*(0:(f_samples-1)).^4+30), ...
    logspace(log10(150),log10(8000),f_samples), ...
    Frequencies));

Loudspeaker_Weights = zeros(f_samples, loudspeakers);
Bright_Zone_sample  = zeros(f_samples, 1);
Quiet_Zone_sample   = zeros(f_samples, 1);


%% Secondly, build a flat spectra desired multizone soundfield for all frequencies from the previous fft and save the speaker signals for each frequency bin.
fprintf('\n====== Building Broadband Multizone Soundfield Speaker Signals ======\n');
fprintf('\tCompletion: ');n=0;
f_rearranged = [floor(f_samples/2):-1:1; floor(f_samples/2)+1:f_samples-rem(f_samples,2)];
f_rearranged = fliplr([f_rearranged(:)' f_samples]);
parfor_progress(f_samples);
parfor f_ = 1:f_samples
    %f_ = f_rearranged( f );
    
    quiet  = Orthogonal_Basis_Expansion.spatial_zone(Frequencies( f_ ), phase, radius, 'quiet');
    bright = Orthogonal_Basis_Expansion.spatial_zone(Frequencies( f_ ), phase, radius, 'pw', 1.0, angle_pw);
    quiet.res  = resolution;
    bright.res = resolution;
    quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');
    bright = bright.setDesiredSoundfield(true, 'suppress_output');
    
    soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
    soundfield = soundfield.addSpatialZone(quiet,  radius2origin, angle2origin(1));
    soundfield = soundfield.addSpatialZone(bright, radius2origin, angle2origin(2));
    
    soundfield.BrightZ_Weight     = 1.0;
    soundfield.QuietZ_Weight      = weight;
    soundfield.UnattendedZ_Weight = 0.5;
    
    soundfield = soundfield.setN( Nplanewaves( f_ ) );
    soundfield = soundfield.createSoundfield('DEBUG', Reproduction_Radius);
    
    
    % Then find the speaker signals required for this multizone soundfield
    setup = Speaker_Setup.loudspeaker_setup;
    setup = setup.addMultizone_Soundfield(soundfield);
    setup.Loudspeaker_Count = loudspeakers;
    setup.Speaker_Arc_Angle = speaker_arc;
    setup.Angle_FirstSpeaker = first_speaker;
    setup = setup.setRadius(speaker_radius);
    
    setup = setup.calc_Loudspeaker_Weights();
    setup = setup.reproduceSoundfield('SAMPLES_ONLY');
    
    
%     Loudspeaker_Weights( f_, :) = setup.Loudspeaker_Weights;
    Bright_Zone_sample( f_ )    = setup.Bright_Sample;%mean(abs( soundfield.Bright_Field(:) ));
    Quiet_Zone_sample( f_ )     = setup.Quiet_Sample;%soundfield.MSE_Quiet_Field ^ (1/2);
    
    %     if ~mod(f,2)
    %         tElapsed = toc;
    %         ratio = f / f_samples;
    %         tRem = (1-ratio) / ratio * tElapsed;
    %         tTot = tElapsed + tRem;
    %         fprintf(repmat('\b',1,n));
    %         n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
    %     end
    parfor_progress;
end
parfor_progress(0);
if ~exist(['+Soundfield_Database\' Output_file_path_ext],'dir'); mkdir(['+Soundfield_Database\' Output_file_path_ext]); end
save(['+Soundfield_Database\' Output_file_path_ext 'Weights_and_Samples__' SetupInfo '.mat'], ...
    'Loudspeaker_Weights', ...
    'Bright_Zone_sample', ...
    'Quiet_Zone_sample' );


%% Finally, apply the speaker weight and reconstruct the loudspeaker signals for each frame of the input signal
load(['+Soundfield_Database\' Output_file_path_ext 'Weights_and_Samples__' SetupInfo '.mat']);

% Here we want to build the speaker signals for each speaker so that our loudspeaker weights are taken into account.
% We want to form the entire spectrum by adding the conjugate of the frame
% to the existing frame where the negative frequencies of the transform
% would usually exist.
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
clear Loudspeaker_Weights; % Save on memory

% We then want to perform an Inverse FFT (ifft) on each full spectrum frame
for frame = 1:size(Loudspeakers_, 1)
%     for spkr = 1:loudspeakers
%         Loudspeakers_(frame,:,spkr) = ifft( Loudspeakers_(frame,:,spkr) );
%     end
    Original(frame,:) = ifft( Original(frame, :) );
end

% Then we should perform the overlap-add method to obtain the complete time domain signal for each speaker
%Loudspeaker_Signals = zeros( [(size(Z,1)+ceil(overlap))*size(Z,2)*2*(1-overlap) loudspeakers] );
% for spkr = 1:loudspeakers
%     Loudspeaker_Signals(:,spkr) = +Tools.OverlapAdd( Loudspeakers_(:,:,spkr), overlap );
% end
Original_ = Tools.OverlapAdd( Original, overlap );
clear Loudspeakers_; % Save on memory

% Normalise Loudspeaker Signals
% Loudspeaker_Signals = Loudspeaker_Signals ./ max(abs(Loudspeaker_Signals(:)));
Original_ = Original_ ./ max(abs(Original_(:)));


% Once we have the speaker signals we should save them for later use as .wav files
% filenumbers = num2str((1:loudspeakers)');
% filenumbers(filenumbers==' ') = '_';
% fullpath = [repmat([Output_file_path Output_file_path_ext Output_file_name], [loudspeakers 1]) ...
%     filenumbers ...
%     repmat(Output_file_ext, [loudspeakers 1]) ];
if ~exist([Output_file_path Output_file_path_ext],'dir'); mkdir([Output_file_path Output_file_path_ext]); end

% for spkr = 1:loudspeakers
%     audiowrite(fullpath(spkr,:), Loudspeaker_Signals(:, spkr), Fs);
% end
audiowrite([Output_file_path Output_file_path_ext ...
    Input_file_name '__Original_' ...
    num2str(f_low ) 'Hz-' ...
    num2str(f_high) 'Hz' ...
    Input_file_ext], Original_, Fs);




%%
Bright_Zone_sample = [zeros(trunc_index_low-1, 1); Bright_Zone_sample; zeros( size(Z,2) - trunc_index_high, 1) ];
Bright_Zone_sample = permute( repmat(Bright_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);

Quiet_Zone_sample  = [zeros(trunc_index_low-1, 1); Quiet_Zone_sample; zeros( size(Z,2) - trunc_index_high, 1) ];
Quiet_Zone_sample  = permute( repmat(Quiet_Zone_sample, [1 1 size(Z,1)]), [3 1 2]);

Bright = zeros( [size(Z,1) (size(Z,2))*2] );
Quiet  = zeros( [size(Z,1) (size(Z,2))*2] );

Bright(:,1:end/2)     = squeeze( Z(:,:,1) ) .* Bright_Zone_sample;
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

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_Results.csv'],'a');
PESQ_Bright = pesq(Original_path,Bright_path);
PESQ_Quiet = pesq(Original_path,Quiet_path);
fprintf(fileID,'%s,PESQ_Bright,%f,PESQ_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_Bright, PESQ_Quiet);
fclose(fileID);

fileID = fopen(['+Results\' Output_file_path_ext 'PESQ_MOS_Results.csv'],'a');
PESQ_MOS_Bright = pesq2mos(PESQ_Bright);
PESQ_MOS_Quiet = pesq2mos(PESQ_Quiet);
fprintf(fileID,'%s,PESQ_MOS_Bright,%f,PESQ_MOS_Quiet,%f,\r\n',[Input_file_name SetupInfo], PESQ_MOS_Bright, PESQ_MOS_Quiet);
fclose(fileID);

%% Calculate and save Active Speech Level values to the results folder
if ~exist(['+Results\' Output_file_path_ext],'dir'); mkdir(['+Results\' Output_file_path_ext]); end

fileID = fopen(['+Results\' Output_file_path_ext 'ASL_Results.csv'],'a');

[ASL_Bright_dB, Activity_Bright] = activlev(Bright_Signal, Fs, 'd');
[ASL_Quiet_dB,  Activity_Quiet ] = activlev(Quiet_Signal,  Fs, 'd');

fprintf(fileID,'%s,ASL_Bright_dB,%f,ASL_Quiet_dB,%f,Activity_Bright,%f,Activity_Quiet,%f,\r\n', ...
        [Input_file_name SetupInfo], ...
        ASL_Bright_dB, ...
        ASL_Quiet_dB, ...
        Activity_Bright, ...
        Activity_Quiet);
fclose(fileID);



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

end