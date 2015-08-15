clc;
clear;
close all;
tic; %Start timing this script
if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))







%% Settings

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
                            % You can suppress the output when executing a complete soundfield construction with 'suppress_output'
Resolution = 150;           % The higher the resolution the better the sampling rate and resultant image
NumberOfPlanewaves = 80;   % The more planewaves the more accurate the sound field (number of linear combined planewaves from 0 degrees to 360 degrees)

                     
Frequencies          = [     150  200  300  400  600  800 ...
                        1000 1500 2000 3000 4000 6000 8000];%[250 500 1000 2000 4000 8000];
Weights              = [0.1 0.2 0.3 0.4 0.6 0.8  ...
                          1   2   3   4   6   8  ...
                          10  20  30  40  60  80 ];
                     
Angles               = 0:90:180;

Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 1.0; % Metres                    

frequencies = length(Frequencies);

%% Results
Contrast__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));
Error_Bright__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Error_Quiet__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));






%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

w=9; %2=10degrees
fprintf('\n====== MSE for Angle Vs Frequency ======\n\n');
fprintf('\tCompletion: ');n=0;
for a = 1:length(Angles)
    parfor f = 1:length(Frequencies)
        %%
        quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies( f ), Radius_of_zones(2), 'quiet');
        bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies( f ), Radius_of_zones(1),   'pw', 1.0, Angles(a));
        quiet.res  = Resolution;
        bright.res = Resolution;
        quiet      = quiet.setDesiredSoundfield(true, 'suppress_output');
        bright     = bright.setDesiredSoundfield(true, 'suppress_output');
        
        %%
        sf = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
        sf = sf.addSpatialZone(quiet,  Position_of_zones(2,1), Position_of_zones(2,2));
        sf = sf.addSpatialZone(bright, Position_of_zones(1,1),   Position_of_zones(1,2));
        
        %%
        sf.QuietZ_Weight = Weights( w );
        sf = sf.setN(NumberOfPlanewaves);
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        
        Contrast__Angle_Vs_Frequency( a, f ) = sf.Acoustic_Brightness_Contrast_dB;
        Error_Bright__Angle_Vs_Frequency( a, f ) = sf.Err_dB_Bright_Field;
        Error_Quiet__Angle_Vs_Frequency( a, f ) = sf.Err_dB_Quiet_Field;
    end
    tElapsed = toc;
    tRem = (1-a/length(Angles)) / (a/length(Angles)) * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', a / length(Angles) * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
end





%%
Title1='Acoustic Brightness Contrast for various Angles & Frequencies';
figure('Name',Title1,'NumberTitle','off')
surf(Frequencies, Angles, Contrast__Angle_Vs_Frequency);
    title(Title1);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Bright Zone Planewave Angle (°)');
    zlabel('Interference Level in Quiet Zone (dB)');
    zlim([-20 120]);
    ZTick = -10:10:120;
    set(gca, 'ZTick', ZTick);

Title2='Loudness in Bright Zone for various Angles & Frequencies';
figure('Name',Title2,'NumberTitle','off')
surf(Frequencies, Angles, Error_Bright__Angle_Vs_Frequency);
    title(Title2);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Bright Zone Planewave Angle (°)');
    zlabel('Interference Level in Quiet Zone (dB)');
    zlim([-120 10]);
    ZTick = -120:10:10;
    set(gca, 'ZTick', ZTick);

Title3='Loudness in Quiet Zone for various Angles & Frequencies';
figure('Name',Title3,'NumberTitle','off')
surf(Frequencies, Angles, Error_Quiet__Angle_Vs_Frequency);
    title(Title3);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Bright Zone Planewave Angle (°)');
    zlabel('Interference Level in Quiet Zone (dB)');
    zlim([-120 10]);
    ZTick = -120:10:10;
    set(gca, 'ZTick', ZTick);





%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script