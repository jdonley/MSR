clc;
%clear;
close all;
tic;
%% Settings

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
Resolution = 100;           % The higher the resolution the better the sampling rate and resultant image

%samples = 120;
samples = 512;
[NumberOfPlanewaves, Frequencies_in_zones] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( ...
                                     samples, ...
                                     28, ...
                                     300, ...
                                     150, ...
                                     8000);
loudness = 20;
SPL = Perceptual_Tools.Loudness( Frequencies_in_zones, loudness );   % SPL of signal in Bright Zone fitting equal loudness curve

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

LUT_resolution = '512f_256w';
angle_pw = 15;
loudspeakers   = 65;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres

%% Generate Quiet Zone Weights so that the results will match the threshold in quiet

% Use the Soundfield_Database Look-Up Table (LUT) for the current
% configuration to get the correct weights
%load +Soundfield_Database\+3m_SpeakerDiameter\Error_Weight_vs_Frequency_15deg.mat
load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat']);

% We need to find the attenuation that is the closest match to the threshold
% The weighting value of 0.05 in the LUT is our reference point (we take
% this to be zero). Because this is an attenuation we should shift all
% attenuation levels for a given frequency so that the 0.05 value gives
% zero dB change.

%Look-Up Table should be for  SPL vs Weights
LUT = Quiet_SPL__Weight_Vs_Frequency - 94;%... %For each frequency set

% Find the weights that provide attenuation just below the threshold
Desired_SPL = Perceptual_Tools.Loudness( Frequencies_in_zones , 0);

% Because the frequencies in the LUT don't match up we will interpolate
% We will then linearly interpolate the closest index that will gives us a weight that
% provides attenuation at the threshold
SPL_Difference = SPL - Desired_SPL;
Frequency_weights = Tools.interpFromVal_2D(LUT, Frequencies, Weights, Frequencies_in_zones, -SPL_Difference);


%% Calculation & Production of Soundfield
fprintf('\n====== Threshold In Quiet Perceptual Weighting - Results ======\n');
fprintf('\tCompletion: ');n=0;
f_rearranged = [(frequencies/2):-1:1; (frequencies/2+1):frequencies];
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
        sf.QuietZ_Weight = Frequency_weights( f, z );
        sf = sf.setN(NumberOfPlanewaves( f ));
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        
        Bright_Error(f) = sf.Err_dB_Bright_Field;
         Quiet_Error(f) = sf.Err_dB_Quiet_Field;    
         
         tElapsed = toc;
         ratio = f_ / frequencies;
         tRem = (1-ratio) / ratio * tElapsed;
         tTot = tElapsed + tRem;
         fprintf(repmat('\b',1,n));
         n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));

    end
end
%%
Bright_SPL = Perceptual_Tools.Loudness( Frequencies, loudness );
Quiet_SPL  = 20 * log10( db2mag(Quiet_Error) .^ (1/2) .* 10.^(SPL' / 20));


%% Results
%load +Perceptual_Models\Error_Limits_120samples.mat
h=figure(1);
subplot(2,2,1)
hold on
plot(Frequencies_in_zones(1, :), Bright_Error);
plot(Frequencies, max(Error_Bright__Weight_Vs_Frequency), ':r', 'LineWidth', 2.0);
plot(Frequencies, min(Error_Bright__Weight_Vs_Frequency), ':', 'Color', [0 2/3 0],'LineWidth', 2.0);
hold off
title('Bright Error', 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
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
plot(Frequencies, Bright_SPL );
title('Bright Sound Pressure Level (SPL)', 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('SPL (dB)','FontSize',20,'FontName','Arial');
ylim([-10 100]);

subplot(2,2,4)
hold on
plot(Frequencies_in_zones(1, :), Quiet_SPL);
plot(Frequencies_in_zones(1, :), Desired_SPL, ':k', 'LineWidth', 2.0);
plot(Frequencies, max(Quiet_SPL__Weight_Vs_Frequency) - 94 + Bright_SPL, ':r', 'LineWidth', 2.0);
plot(Frequencies, min(Quiet_SPL__Weight_Vs_Frequency) - 94 + Bright_SPL, ':', 'Color', [0 2/3 0],'LineWidth', 2.0);
hold off
title('Quiet Sound Pressure Level (SPL)', 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('SPL (dB)','FontSize',20,'FontName','Arial');
ylim([-10 100]);

%%
%saveas(h,'+Perceptual_Models\Constraint_on_Quiet__ThresholdInQuiet.fig');


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
