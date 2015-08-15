clc;
%clear;
close all;
tic;
%% Settings

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
Resolution = 100;           % The higher the resolution the better the sampling rate and resultant image

samples = 120;
Frequencies_in_zones = logspace(log10(150),log10(8000),samples);   % Zone frequencies

quiet_mask_noise_dB = 45; % Phons (dB)

loudness = 60;
SPL = Perceptual_Tools.Loudness( Frequencies_in_zones, loudness );   % SPL of signal in Bright Zone fitting equal loudness curve
%SPL = ones(length(Frequencies_in_zones),1) * loudness; % Flat response

NumberOfPlanewaves   = ceil((300-30)/((samples-1)^4)*(0:(samples-1)).^4+30);   % The more planewaves the more accurate the sound field


Angles_in_zones      = [15; ...         % Planewave angle Zone 1
                        90];           % Planewave angle Zone 2
Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 1.0; % Metres                    

[zones, frequencies] = size(Frequencies_in_zones);  % Currently only support for two zones using Orthogonal Basis Expansion


sf = [];
Soundfield = [];
BrightField = [];

Bright_Error = zeros(length(Frequencies_in_zones),1);
Quiet_Error  = zeros(length(Frequencies_in_zones),1);


%% Generate Quiet Zone Weights so that the results will match the threshold in quiet

% Use the Soundfield_Database Look-Up Table (LUT) for the current
% configuration to get the correct weights
load +Soundfield_Database\+3m_SpeakerDiameter\Error_Weight_vs_Frequency_15deg.mat

% We need to find the attenuation that is the closest match to the threshold
% The weighting value of 0.05 in the LUT is our reference point (we take
% this to be zero). Because this is an attenuation we should shift all
% attenuation levels for a given frequency so that the 0.05 value gives
% zero dB change.

%Look-Up Table should be for  SPL vs Weights
LUT = Quiet_SPL__Weight_Vs_Frequency(: , 1:end) - 94;%... %For each frequency set

% Find the weights that provide attenuation just below the threshold
%Desired_SPL = Perceptual_Tools.Threshold_in_Quiet( Frequencies_in_zones(1, :) , 'ISO226'); % ISO226 Standard
Desired_SPL = Perceptual_Tools.Loudness( Frequencies_in_zones(1, :) , quiet_mask_noise_dB); 

% Because the frequencies in the LUT don't match up we will interpolate
% We will linearly interpolate the closest index that will gives us a weight that
% provides attenuation at the threshold
SPL_Difference = Desired_SPL - SPL;
Frequency_weights(1:zones, :) = Tools.interpFromVal_2D(LUT, Frequencies, Weights, Frequencies_in_zones(1,:), SPL_Difference);

%% Calculation & Production of Soundfield
fprintf('\n====== Threshold In Quiet Perceptual Weighting - Results ======\n');
fprintf('\tCompletion: ');n=0;
f_rearranged = [(frequencies/4):-1:1; (frequencies*3/4+1):frequencies; (frequencies/4+1):frequencies/2; frequencies*3/4:-1:(frequencies/2+1)];
f_rearranged = fliplr(f_rearranged(:)');
for z = 1:zones
    for f_ = 1:frequencies
        f = f_rearranged(f_);
        %%
        quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones( f ), 0, Radius_of_zones(3-z), 'quiet');
        bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones( f ), 0, Radius_of_zones(z),   'pw', 1.0, Angles_in_zones(z));
        quiet.res  = Resolution;
        bright.res = Resolution;
        quiet      =  quiet.setDesiredSoundfield(true, 'suppress_output');
        bright     = bright.setDesiredSoundfield(true, 'suppress_output');

        
        %%
        sf = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
        sf = sf.addSpatialZone(quiet,  Position_of_zones(3-z,1), Position_of_zones(3-z,2));
        sf = sf.addSpatialZone(bright, Position_of_zones(z,1),   Position_of_zones(z,2));
        
        %%
        sf.QuietZ_Weight = Frequency_weights(z, f);
        sf = sf.setN(NumberOfPlanewaves( f ));
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        sf = sf.norm_soundfield;        
        
        Bright_Error(f) = sf.Err_dB_Bright_Field;
        Quiet_Error(f) = sf.Err_dB_Quiet_Field;
        if ~mod(f_,4)
            tElapsed = toc;
            ratio = f_ / frequencies;
            tRem = (1-ratio) / ratio * tElapsed;
            tTot = tElapsed + tRem;
            fprintf(repmat('\b',1,n));
            n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
        end
    end
end


%% Results
Bright_SPL = SPL;
Quiet_SPL  = 20 * log10( db2mag(Quiet_Error) .^ (1/2) .* 10.^(SPL / 20));

load +Perceptual_Models\Error_Limits_120samples.mat

Bright_Improvement = mean(Bright_Max_Error - Bright_Error);
Quiet_Perceivable  = mean(Quiet_SPL - Desired_SPL);

%% Plots
h=figure(1);
subplot(2,2,1)
hold on
plot(Frequencies_in_zones(1, :), Bright_Error);
plot(Error_Limits_Frequencies, Bright_Max_Error, ':r', 'LineWidth', 2.0);
plot(Error_Limits_Frequencies, Bright_Min_Error, ':', 'Color', [0 2/3 0],'LineWidth', 2.0);
hold off
title(['Bright Error' 10 '(' num2str(round(Bright_Improvement*100)/100) 'dB less error)'], 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('Error (dB)','FontSize',20,'FontName','Arial');
ylim([-80 0]);

subplot(2,2,2)
plot(Frequencies_in_zones(1, :), Quiet_Error);
title('Quiet Error', 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('Error (dB)','FontSize',20,'FontName','Arial');
ylim([-80 0]);

subplot(2,2,3)
plot(Frequencies_in_zones(1, :), Bright_SPL );
title(['Bright SPL (Loudness: ' num2str(loudness) 'dB)'], 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('SPL (dB)','FontSize',20,'FontName','Arial');
ylim([-10 100]);

subplot(2,2,4)
hold on
plot(Frequencies_in_zones(1, :), Quiet_SPL);
plot(Frequencies_in_zones(1, :), Desired_SPL, ':k', 'LineWidth', 2.0);
plot(Frequencies_in_zones(1, :), Perceptual_Tools.Threshold_in_Quiet( Frequencies_in_zones, 'ISO226' ), '-k', 'LineWidth', 1.0);
plot(Error_Limits_Frequencies, Quiet_Max_SPL+(Bright_SPL-Perceptual_Tools.Loudness( Frequencies_in_zones, 20 )), ':r', 'LineWidth', 2.0);
plot(Error_Limits_Frequencies, Quiet_Min_SPL+(Bright_SPL-Perceptual_Tools.Loudness( Frequencies_in_zones, 20 )), ':', 'Color', [0 2/3 0],'LineWidth', 2.0);
hold off
title(['Quiet SPL (Mask Noise: ' num2str(quiet_mask_noise_dB) 'dB)' 10 '(Perceivable: ' num2str(round(Quiet_Perceivable*100)/100) 'dB)'], 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('SPL (dB)','FontSize',20,'FontName','Arial');
ylim([-10 100]);

%%
saveas(h,['+Perceptual_Models\Constraint_on_Quiet__ISO226_ComfortNoise_' num2str(quiet_mask_noise_dB) 'dB.fig']);


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
