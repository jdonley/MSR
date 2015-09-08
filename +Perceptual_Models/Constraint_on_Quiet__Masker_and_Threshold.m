clc;
%clear;
close all;
tic;
%% Settings

%Figure Output Settings
DocumentPath = 'tex\latex\APSIPA';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10 + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 5/3;
FontSize = 9;
Font = 'Times';

% Simulation resolution samples per metre
DebugMode            = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
Resolution           = 100;           % The higher the resolution the better the sampling rate and resultant image


samples              = 512 ;
[NumberOfPlanewaves, Frequencies_in_zones] = Soundfield_Database.LUT_Builders.Orthogonal_Planewave_Selection( ...
                                     samples, ...
                                     28, ...
                                     300, ...
                                     150, ...
                                     8000);

loudness = 30;
SPL = Perceptual_Tools.Loudness( Frequencies_in_zones, loudness )';   % SPL of signal in Bright Zone fitting equal loudness curve
%SPL = ones(length(Frequencies_in_zones),1) * loudness; % Flat response

Masker_frequency     = 2000; % Hz
Masker_level         = loudness;  % dB


Angles_in_zones      = [15; ...         % Planewave angle Zone 1
                        15];           % Planewave angle Zone 2
Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius  = 1.0; % Metres

frequencies          = length(Frequencies_in_zones);  % Currently only support for two zones using Orthogonal Basis Expansion

zones = 1;

sf = [];
Soundfield = [];
BrightField = [];

Bright_Error = zeros(length(Frequencies_in_zones),1);
Quiet_Error  = zeros(length(Frequencies_in_zones),1);
Quiet_SPL    = zeros(length(Frequencies_in_zones),1);

LUT_resolution = '512f_256w';
angle_pw = 15;
loudspeakers   = 65;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres

%% Generate Quiet Zone Weights so that the results will match the threshold in quiet

% Use the Soundfield_Database Look-Up Table (LUT) for the current
% configuration to get the correct weights
load(['+Soundfield_Database\+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc\LUT_Weight_vs_Frequency_' num2str(angle_pw) 'deg_' LUT_resolution '.mat']);

% We need to find the attenuation that is the closest match to the threshold
% The weighting values were built from a desried SPL of 94dB as reference,
% we then attenuate the quiet zone at each frequency to match the desired SPL.


%Look-Up Table should be for  SPL vs Weights
LUT = Quiet_SPL__Weight_Vs_Frequency(: , 1:end) - 94;%... %For each frequency set


% Find the weights that provide attenuation just below the combined
% threshold and masker.
Thresh_  = Perceptual_Tools.Threshold_in_Quiet( Frequencies_in_zones , 'ISO226' )'; % Accepts frequency in kHz
z_maskee = Perceptual_Tools.FrequencyToBark( Frequencies_in_zones );
z_masker = Perceptual_Tools.FrequencyToBark( Masker_frequency );
Mask_    = Perceptual_Tools.Masking_Curves.ISO_IEC_MPEG_Model2( z_maskee, z_masker, Masker_level );
Desired_SPL = Perceptual_Tools.Combine_ThreshMask( Thresh_, Mask_ );


% Because the frequencies in the LUT don't match up we will interpolate
% We will then linearly interpolate the closest index that will gives us a weight that
% provides attenuation at the threshold
SPL_Difference = SPL - Desired_SPL;
Frequency_weights = Tools.interpFromVal_2D(LUT, Frequencies, Weights, Frequencies_in_zones, -SPL_Difference);


%% Calculation & Production of Soundfield
fprintf('\n====== %.2fkHz Masker & Threshold Perceptual Weighting - Results ======\n', Masker_frequency/1000);
fprintf('\tCompletion: ');n=0;
f_rearranged = [(frequencies/2):-1:1; (frequencies/2+1):frequencies];
f_rearranged = fliplr(f_rearranged(:)');
z = zones;
for f_ = 1:frequencies
    f = f_rearranged(f_);
    %%
    quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones( f ), 0, Radius_of_zones( 3-z ), 'quiet');
    bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones( f ), 0, Radius_of_zones( z ),   'pw', 1.0, Angles_in_zones( z ));
    quiet.res  = Resolution;
    bright.res = Resolution;
    quiet      =  quiet.setDesiredSoundfield(true, 'suppress_output');
    bright     = bright.setDesiredSoundfield(true, 'suppress_output');
    
    %% If we come across a masker frequency at the same frequency as in the bright zone then we do an improved synthesis
     % We replace the desired quite zone with the masker frequency and
     % match it's phase with the bright zone
     if Masker_frequency == Frequencies_in_zones( f )
         perpendicular_dist = Radius_of_zones( z ) + Radius_of_zones( 3-z ); % TODO make this the true perpendicular distance between zones
         wavelength = 343 / Frequencies_in_zones( f );
         phase = -rem( perpendicular_dist * cos(Angles_in_zones( z )), wavelength) / wavelength * 360;
         temp = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones( f ), phase, Radius_of_zones( 3-z ),'pw', 1.0, Angles_in_zones( z ));
         temp.res = Resolution;
         temp = temp.setDesiredSoundfield(true, 'suppress_output');
         quiet.Soundfield_d = temp.Soundfield_d;
     end
    %%
    sf = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
    sf = sf.addSpatialZone(quiet,  Position_of_zones(3-z,1), Position_of_zones(3-z,2));
    sf = sf.addSpatialZone(bright, Position_of_zones(z,1),   Position_of_zones(z,2));
    
    %%
    sf.QuietZ_Weight = Frequency_weights( f, z );
    sf = sf.setN(NumberOfPlanewaves(z, f));
    sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
    
    Bright_Error(f) = sf.Err_dB_Bright_Field;
    Quiet_Error(f)  = sf.Err_dB_Quiet_Field;
    Quiet_SPL(f)  = 20 * log10( db2mag(Quiet_Error(f)) ^ (1/2) * 10^(SPL(f) / 20));

    % If we have adjusted because of the masker frequency then we need to
    % recalculate the error as if we have more than one bright zone
     if Masker_frequency == Frequencies_in_zones( f )
         sf.Bright_Zone.Soundfield_d = sf.Quiet_Zone.Soundfield_d;
         sf.Bright_Zone.Soundfield_d_mask = sf.Quiet_Zone.Soundfield_d_mask;
         sf.Bright_Field = sf.Quiet_Field;
         sf.ErrorInBright_Calc();
         Quiet_Error(f) = sf.Err_dB_Bright_Field;
         Quiet_SPL(f) = Masker_level;
     end    
         
    tElapsed = toc;
    ratio = f_ / frequencies;
    tRem = (1-ratio) / ratio * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
    
end
%%
Bright_SPL = Perceptual_Tools.Loudness( Frequencies, loudness )';


%% Results
close all;
h=figure(1);
N_axes = 3;

subplot(1,N_axes,1)
t(1) = title('Bright Zone SPL', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
hold on;
set(gca, 'ColorOrderIndex',6); pl(1) = plot(Frequencies, Bright_SPL , 'LineWidth', 1.0);
pl(2) = plot(Frequencies, Perceptual_Tools.Loudness( Frequencies, 0 ), '--k', 'LineWidth', 0.5);
hold off;
xlabel('Frequency (Hz)');
ylabel('SPL (dB)');
ylim([-10 60]);
grid on;box on;%grid minor;grid minor;
set(gca,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1], ...
            'XScale', 'log');
le = legend(pl(1:2),{'Desired SPL';'Threshold in Quiet'},'Location','northeast','FontSize',ceil(FontSize-2),'Box','off');
le.Position = le.Position + [0.028 0.19 0 0];

subplot(1,N_axes,2)
t(2) = title('Bright Zone Spatial Error', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
hold on
set(gca, 'ColorOrderIndex',6); pl(3) = plot(Frequencies_in_zones(1, :), Bright_Error, 'LineWidth', 1.0);
set(gca, 'ColorOrderIndex',7); pl(4) = plot(Frequencies, max(Error_Bright__Weight_Vs_Frequency), '--', 'LineWidth', 0.5);
set(gca, 'ColorOrderIndex',5); pl(5) = plot(Frequencies, min(Error_Bright__Weight_Vs_Frequency), '--','LineWidth', 0.5);
hold off
xlabel('Frequency (Hz)');
ylabel('Error (dB)');
ylim([-60 10]);
grid on;box on;%grid minor;grid minor;
set(gca,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1], ...
            'XScale', 'log');
le = legend(pl([4 3 5]),{'Max Error';'Actual Error';'Min Error'},'Location','southwest','FontSize',ceil(FontSize-2),'Box','off');
le.Position = le.Position + [-0.050 -0.205 0 0];

subplot(1,N_axes,3)
t(3) = title('Quiet Zone SPL', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
hold on
set(gca, 'ColorOrderIndex',6); pl(6) = plot(Frequencies_in_zones(1, :), Quiet_SPL, 'LineWidth', 1.0);
pl(7) = plot(Frequencies_in_zones(1, :), Desired_SPL, '-.k', 'LineWidth', 0.5);
set(gca, 'ColorOrderIndex',7); pl(8) = plot(Frequencies, max(Quiet_SPL__Weight_Vs_Frequency) - 94 + Bright_SPL, '--', 'LineWidth', 0.5);
set(gca, 'ColorOrderIndex',5); pl(9) = plot(Frequencies, min(Quiet_SPL__Weight_Vs_Frequency) - 94 + Bright_SPL, '--', 'LineWidth', 0.5);
hold off
xlabel('Frequency (Hz)');
ylabel('SPL (dB)');
ylim([-10 60]);
grid on;box on;%grid minor;grid minor;
set(gca,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1], ...
            'XScale', 'log');
le = legend(pl([8 6 7 9]),{'Max SPL';'Leaked SPL';'Desired SPL';'Min SPL'},'Location','northeast','FontSize',ceil(FontSize-2),'Box','off');
le.Position = le.Position + [0.118 0.156 0 0];
        
tightfig;
factor1 = 2.2;
factor2 = 0.25;
factor3 = -0.4;
title_offset = FontSize/4;
offset_cm = FontSize*254/720 /10 * factor1;
for i=1:N_axes, t(i).Position(2) = t(i).Position(2)+title_offset;end
set(gcf, 'PaperUnits','centimeters', 'PaperSize', [plot_width+factor3, plot_width/aspect_ratio/N_axes + offset_cm],'PaperPosition',[factor3 factor2 plot_width,plot_width/aspect_ratio/N_axes + offset_cm]);
%tightfig;
%print(['-d' print_fmt], ['-r' num2str(print_res)], ['out.' print_fmt]);
print(['-d' print_fmt], [DocumentPath '\Psychoacoustic_Error_Reduction']);

% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;
%%
%saveas(h,'+Perceptual_Models\Constraint_on_Quiet__3kHzMaskerThreshold.fig');


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
