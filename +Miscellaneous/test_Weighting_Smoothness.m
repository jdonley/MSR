clc;
%clear;
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
Resolution = 100;           % The higher the resolution the better the sampling rate and resultant image
% Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.

SPL = 94; % 94dB is 1Pa (RMS) when reference to 20uPa (RMS)

Number_of_Frequencies = 1;
Number_of_Weights = 28;

NumberOfPlanewaves   =       50;   % The more planewaves the more accurate the sound field (more planewaves increases error from matrix inversion, less planewaves increases spatial aliasing error)
Frequencies          =      500 ;

Weights              =      logspace(log10(0.05), log10(  90), Number_of_Weights    ) ;
                     
Angles               = [15    ];

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
Quiet_SPL__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));


%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

fprintf('\n====== MSE for Frequency Vs Weight ======\n\n');
fprintf('\tCompletion: ');n=0;
a=1;
for w = 1:length(Weights)
    parfor f = 1:length(Frequencies)
        %%
        quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies( f ), 0, Radius_of_zones(2), 'quiet');
        bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies( f ), 0, Radius_of_zones(1),   'pw', 1.0, Angles(a));
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
        Quiet_SPL__Weight_Vs_Frequency( w, f )  = 20 * log10( sf.MSE_Quiet_Field .^ (1/2) * 10^(SPL / 20));
    end
    tElapsed = toc;
    ratio = w/length(Weights);
    tRem = (1-ratio) / ratio * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
end





%%
Title1='Acoustic Brightness Contrast for various Weights & Frequencies';
h=figure('Name',Title1,'NumberTitle','off');
plot(Weights, Contrast__Weight_Vs_Frequency);
    title(Title1);
    xlabel('Quiet Zone Weight');
    set(gca, 'XScale', 'log');
    ylabel('Acoustic Brightness Contrast (dB)');
    ylim([-20 120]);
    YTick = -10:10:120;
    set(gca, 'YTick', YTick);
    
Title2='Error in Bright Zone for various Weights & Frequencies';
h=figure('Name',Title2,'NumberTitle','off');
plot(Weights, Error_Bright__Weight_Vs_Frequency);
    title(Title2);
    xlabel('Quiet Zone Weight');
    set(gca, 'XScale', 'log');
    ylabel('Mean Squared Error in Bright Zone (dB)');
    ylim([-120 10]);
    YTick = -120:10:10;
    set(gca, 'YTick', YTick);

Title3='Error in Quiet Zone for various Weights & Frequencies';
h=figure('Name',Title3,'NumberTitle','off');
plot(Weights, Error_Quiet__Weight_Vs_Frequency);
    title(Title3);
    xlabel('Quiet Zone Weight');
    set(gca, 'XScale', 'log');
    ylabel('Mean Squared Error in Quiet Zone (dB)');
    ylim([-120 10]);
    YTick = -120:10:10;
    set(gca, 'YTick', YTick);
    
Title4='SPL in Quiet Zone for various Weights & Frequencies';
h=figure('Name',Title4,'NumberTitle','off');
plot(Weights, Quiet_SPL__Weight_Vs_Frequency);
    title(Title4);
    xlabel('Quiet Zone Weight');
    set(gca, 'XScale', 'log');
    ylabel('SPL in Quiet Zone (dB)');
    ylim([-10 120]);
    YTick = -10:10:120;
    set(gca, 'YTick', YTick);



%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script