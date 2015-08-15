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




NumberOfPlanewaves   = repmat(10:2:300, [20 length(10:2:300)]) 
                    
Frequencies          = logspace(log10(100), log10(8000), 20);

Weight               = 2.5;
Angle                = 90;

Reproduction_Radius = 1.0; % Metres                    

[frequencies N] = size(NumberOfPlanewaves);

%% Results
Contrast__Frequency_Vs_N     = zeros( size(NumberOfPlanewaves) );
Error_Bright__Frequency_Vs_N = zeros( size(NumberOfPlanewaves) );
Error_Quiet__Frequency_Vs_N  = zeros( size(NumberOfPlanewaves) );






%% Calculation of MSE for Weight Vs Angle (MSE__Weight_Vs_Angle)

fprintf('\n====== MSE for Frequency Vs Number of Linear Combined Planewaves ======\n\n');
fprintf('\tCompletion: ');nt=0;
for n = 1:N
    parfor f = 1:frequencies
        %%
        quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies( f ), 0, 0.3, 'quiet');
        bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies( f ), 0, 0.3,   'pw', 1.0, Angle);
        quiet.res  = Resolution;
        bright.res = Resolution;
        quiet      = quiet.setDesiredSoundfield(true, 'suppress_output');
        bright     = bright.setDesiredSoundfield(true, 'suppress_output');
        
        %%
        sf = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
        sf = sf.addSpatialZone(quiet,  0.6, 0);
        sf = sf.addSpatialZone(bright, 0.6, 180);
        
        %%
        sf.QuietZ_Weight = Weight;
        sf = sf.setN( NumberOfPlanewaves( f, n ) );
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        
        Contrast__Frequency_Vs_N( f, n ) = sf.Acoustic_Brightness_Contrast_dB;
        Error_Bright__Frequency_Vs_N( f, n ) = sf.Err_dB_Bright_Field;
        Error_Quiet__Frequency_Vs_N( f, n ) = sf.Err_dB_Quiet_Field;
    end
    tElapsed = toc;
    tRem = (1-n/N) / (n/N) * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,nt));
    nt=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', n / N * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
end





%%


Title='Loudness in Quiet Zone for various Frequencies & Number of Planewaves';
figure('Name',Title,'NumberTitle','off') 
surf( Error_Quiet__Frequency_Vs_N' );
legend show Location NorthEastOutside;
legend(['100Hz';'150Hz';'200Hz';'300Hz';'400Hz';'600Hz';'800Hz';'1kHz ';]);

save Error_Quiet__100hz_1khz_Vs_N.mat Error_Quiet__Frequency_Vs_N


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
