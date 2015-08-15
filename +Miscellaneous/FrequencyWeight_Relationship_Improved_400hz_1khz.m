clc;
%clear;
close all;
tic; %Start timing this script
[ST,I]=dbstack;
disp(ST.name);
if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))







%% Settings

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
                            % You can suppress the output when executing a complete soundfield construction with 'suppress_output'
Resolution = 150;           % The higher the resolution the better the sampling rate and resultant image

NumberOfPlanewaves   = [  20   22   22   26   38   38   48   72   72   94   94 120 ...
                         120  120  ];%120  120  120  120  120  120  120  120  120]; % Minimum values with least error caused by aliasing
Frequencies          = [ 100  125  150  200  250  300  400  500  600  700  800 900 ...
                        1000 1250 ];%1500 2000 2500 3000 4000 5000 6000 7000 8000];
Weights              = [0.1 0.2 0.3 0.4 0.6 0.8  ...
                          1   2   3   4   6   8  ...
                          10  20  30  40  60  80 ];
                     
Angles               = [14    ];

Radius_of_zones      = [0.9; ...        % Zone 1 radius
                        0.9];           % Zone 2 radius
Position_of_zones    = [1.8, 180; ...   % Zone 1 distance and angle from origin
                        1.8, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 3.0; % Metres                    

frequencies = length(Frequencies);

%% Results
Contrast__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));
Error_Bright__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Error_Quiet__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));






%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

fprintf('\n====== MSE for Frequency Vs Weight ======\n\n');
fprintf('\tCompletion: ');n=0;
a=1;
for w = 1:length(Weights)
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
        sf = sf.setN( NumberOfPlanewaves( f ) );
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        
        Contrast__Weight_Vs_Frequency( w, f ) = sf.Acoustic_Brightness_Contrast_dB;
        Error_Bright__Weight_Vs_Frequency( w, f ) = sf.Err_dB_Bright_Field;
        Error_Quiet__Weight_Vs_Frequency( w, f ) = sf.Err_dB_Quiet_Field;
    end
    tElapsed = toc;
    tRem = (1-w/length(Weights)) / (w/length(Weights)) * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', w / length(Weights) * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
end





%%
Title1='Acoustic Brightness Contrast for various Weights & Frequencies';
h=figure('Name',Title1,'NumberTitle','off');
surf(Frequencies, Weights, Contrast__Weight_Vs_Frequency);
    title(Title1);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Quiet Zone Weight');
    set(gca, 'YScale', 'log');
    zlabel('Acoustic Brightness Contrast (dB)');
    zlim([-20 120]);
    ZTick = -10:10:120;
    set(gca, 'ZTick', ZTick);
    saveas(h,[Title1 ' 400hz-1khz.fig']);
    
Title2='Loudness in Bright Zone for various Weights & Frequencies';
h=figure('Name',Title2,'NumberTitle','off');
surf(Frequencies, Weights, Error_Bright__Weight_Vs_Frequency);
    title(Title2);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Quiet Zone Weight');
    set(gca, 'YScale', 'log');
    zlabel('Mean Squared Error in Bright Zone (dB)');
    zlim([-120 10]);
    ZTick = -120:10:10;
    set(gca, 'ZTick', ZTick);
    saveas(h,[Title2 ' 400hz-1khz.fig']);

Title3='Loudness in Quiet Zone for various Weights & Frequencies';
h=figure('Name',Title3,'NumberTitle','off');
surf(Frequencies, Weights, Error_Quiet__Weight_Vs_Frequency);
    title(Title3);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Quiet Zone Weight');
    set(gca, 'YScale', 'log');
    zlabel('Mean Squared Error in Quiet Zone (dB)');
    zlim([-120 10]);
    ZTick = -120:10:10;
    set(gca, 'ZTick', ZTick);
    saveas(h,[Title3 ' 400hz-1khz.fig']);





%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script