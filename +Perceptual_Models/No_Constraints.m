clc;
%clear;
%close all;
tic;
%% Settings

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
Resolution = 100;           % The higher the resolution the better the sampling rate and resultant image


samples              = 120;
Frequencies_in_zones = [ logspace(log10(150),log10(8000),samples);... % Zone 1 frequencies
                         logspace(log10(150),log10(8000),samples) ];   % Zone 2 frequencies
                     
NumberOfPlanewaves   = [ ceil((300-30)/((samples-1)^4)*(0:(samples-1)).^4+30)   ; ...
                         ceil((300-30)/((samples-1)^4)*(0:(samples-1)).^4+30)   ];   % The more planewaves the more accurate the sound field
                     
loudness = 20;
SPL = Perceptual_Tools.Loudness( Frequencies_in_zones(1,:), loudness );   % SPL of signal in Bright Zone fitting equal loudness curve
%SPL = ones(length(Frequencies_in_zones),1) * loudness; % Flat response
                     
                     
Frequency_weights    = [ linspace( 0.05, 0.05, length(Frequencies_in_zones))   ;...     % Importance weight for quiet zone of Zone 1
                         linspace( 0.05, 0.05, length(Frequencies_in_zones))   ];       % Importance weight for quiet zone of Zone 2
                     
Angles_in_zones      = [15; ...         % Planewave angle Zone 1
                        90];           % Planewave angle Zone 2
Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 1.0; % Metres                    

[zones, frequencies] = size(Frequencies_in_zones);  % Currently only support for two zones using Orthogonal Basis Expansion

zones = 1; % Force just one zone to be created. (Comment out for creation of both zones)

sf = [];
Soundfield = [];
BrightField = [];

Bright_Error = zeros(length(Frequencies_in_zones),1);
Quiet_Error  = zeros(length(Frequencies_in_zones),1);


%% Calculation & Production of Soundfield
fprintf('\n====== No Perceptual Weighting - Results ======\n');
fprintf('\tCompletion: ');n=0;
f_rearranged = [(frequencies/2):-1:1; (frequencies/2+1):frequencies];
f_rearranged = fliplr(f_rearranged(:)');
for z = 1:zones
    for f_ = 1:frequencies
        f = f_rearranged(f_);
        %%
        quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones(z, f), 0, Radius_of_zones(3-z), 'quiet');
        bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones(z, f), 0, Radius_of_zones(z),   'pw', 1.0, Angles_in_zones(z));
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
        sf = sf.setN(NumberOfPlanewaves(z, f));
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        sf = sf.norm_soundfield;        
        
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
Bright_SPL = SPL;
Quiet_SPL  = 20 * log10( db2mag(Quiet_Error) .^ (1/2) .* 10.^(SPL / 20));


%% Results
h=figure(1);
subplot(2,2,1)
plot(Frequencies_in_zones(1, :), Bright_Error);
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
plot(Frequencies_in_zones(1, :), Bright_SPL );
title('Bright Sound Pressure Level (SPL)', 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('SPL (dB)','FontSize',20,'FontName','Arial');
ylim([-10 100]);

subplot(2,2,4)
plot(Frequencies_in_zones(1, :), Quiet_SPL);
title('Quiet Sound Pressure Level (SPL)', 'FontWeight','bold', 'FontSize',24, 'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Arial');
set(gca, 'XScale', 'log');
ylabel('SPL (dB)','FontSize',20,'FontName','Arial');
ylim([-10 100]);

%%
Bright_Min_Error = Bright_Error;
Quiet_Min_SPL    = Quiet_SPL;
Error_Limits_Frequencies = Frequencies_in_zones(1, :);
save(['+Perceptual_Models\Error_Limits_' num2str(samples) 'samples.mat'], ...
    'Bright_Min_Error', ...
    'Quiet_Min_SPL', ...
    'Error_Limits_Frequencies', ...
    '-append');
%%
saveas(h,'+Perceptual_Models\No_Constraints.fig');

%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
