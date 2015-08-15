clc;
%clear;
close all;
tic; %Start timing this script
if (matlabpool('size') ~= 0); matlabpool close; end; %Close existing pool if open
matlabpool %Start new pool
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))







%% Settings
Database_Path = '+Soundfield_Database\+3m_SpeakerDiameter\';

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
                            % You can suppress the output when executing a complete soundfield construction with 'suppress_output'
Resolution = 100;           % The higher the resolution the better the sampling rate and resultant image
% Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.

SPL = 94; % 94dB is 1Pa (RMS) when reference to 20uPa (RMS)

Number_of_Frequencies = 64;
Number_of_Weights = 16;

Starting_Frequency  =  150; %Hz
Starting_No_of_Pw   =   30;
Finishing_Frequency = 8000; %Hz
Finishing_No_of_Pw   =  300;

NumberOfPlanewaves   = ceil((Finishing_No_of_Pw-Starting_No_of_Pw)/((Number_of_Frequencies-1)^4)*(0:(Number_of_Frequencies-1)).^4+Starting_No_of_Pw);   % The more planewaves the more accurate the sound field
Frequencies          = logspace(log10(Starting_Frequency), log10(Finishing_Frequency), Number_of_Frequencies) ;

Weights              = [0     logspace(log10( 1e-2), log10(  1e4), Number_of_Weights - 1) ];
                     
Angles               = [15    ];

Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 1.0; % Metres  

% Speaker setup variables
loudspeakers   = 16;  % Number of loudspeakers
speaker_arc    = 180;  % Degrees
first_speaker  = 90+15; % Degrees
speaker_radius = 1.5; % Metres

frequencies = length(Frequencies);

%% Results
Contrast__Weight_Vs_Frequency     = zeros(length(Weights), length(Frequencies));
Error_Bright__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Error_Quiet__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));
Quiet_SPL__Weight_Vs_Frequency    = zeros(length(Weights), length(Frequencies));

%Samples
Bright_Sample__Weight_Vs_Frequency = zeros(length(Weights), length(Frequencies));
Quiet_Sample__Weight_Vs_Frequency  = zeros(length(Weights), length(Frequencies));


%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

fprintf('\n====== Building Look-Up Table (Frequency Weights) ======\n\n');
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
        
        
        setup = Speaker_Setup.loudspeaker_setup;
        setup = setup.addMultizone_Soundfield(sf);
        setup.Loudspeaker_Count = loudspeakers;
        setup.Speaker_Arc_Angle = speaker_arc;
        setup.Angle_FirstSpeaker = first_speaker;
        setup = setup.setRadius(speaker_radius);
    
        setup = setup.calc_Loudspeaker_Weights();
        setup = setup.reproduceSoundfield('SAMPLES_ONLY');
    
    
        %Loudspeaker_Weights( f_, :) = setup.Loudspeaker_Weights;
        Bright_Sample__Weight_Vs_Frequency( w, f ) = setup.Bright_Sample;
        Quiet_Sample__Weight_Vs_Frequency( w, f ) = setup.Quiet_Sample;
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
    saveas(h,[Database_Path Title1 '.fig']);
    
Title2='Error in Bright Zone for various Weights & Frequencies';
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
    saveas(h,[Database_Path Title2 '.fig']);

Title3='Error in Quiet Zone for various Weights & Frequencies';
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
    saveas(h,[Database_Path Title3 '.fig']);
    
Title4='SPL in Quiet Zone for various Weights & Frequencies';
h=figure('Name',Title4,'NumberTitle','off');
surf(Frequencies, Weights, Quiet_SPL__Weight_Vs_Frequency);
    title(Title4);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Quiet Zone Weight');
    set(gca, 'YScale', 'log');
    zlabel('SPL in Quiet Zone (dB)');
    zlim([-10 120]);
    ZTick = -10:10:120;
    set(gca, 'ZTick', ZTick);
    saveas(h,[Database_Path Title4 '.fig']);
    
%% Samples 
Title5='Sample in Bright Zone for various Weights & Frequencies';
h=figure('Name',Title5,'NumberTitle','off');
surf(Frequencies, Weights, abs(Bright_Sample__Weight_Vs_Frequency));
    title(Title5);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Quiet Zone Weight');
    set(gca, 'YScale', 'log');
    zlabel('Bright Zone Sample');
    saveas(h,[Database_Path Title5 '.fig']);
    
Title6='Sample in Quiet Zone for various Weights & Frequencies';
h=figure('Name',Title6,'NumberTitle','off');
surf(Frequencies, Weights, abs(Quiet_Sample__Weight_Vs_Frequency));
    title(Title6);
    xlabel('Frequency (Hz)');
    set(gca, 'XScale', 'log');
    ylabel('Quiet Zone Weight');
    set(gca, 'YScale', 'log');
    zlabel('Bright Zone Sample');
    saveas(h,[Database_Path Title6 '.fig']);

%%
save([Database_Path 'LUT_Weight_vs_Frequency_' ...
      num2str(Angles) 'deg_' ...
      num2str(Number_of_Frequencies) 'f_' ...
      num2str(Number_of_Weights) 'w.mat'], ...
     'Frequencies', ...
     'Weights', ...
     'Error_Quiet__Weight_Vs_Frequency', ...
     'Error_Bright__Weight_Vs_Frequency', ...
     'Contrast__Weight_Vs_Frequency', ...
     'Quiet_SPL__Weight_Vs_Frequency', ...
     'Bright_Sample__Weight_Vs_Frequency', ...
     'Quiet_Sample__Weight_Vs_Frequency');


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script