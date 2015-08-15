clc;
%clear;
%close all;
tic;
%% Settings

% Simulation resolution samples per metre
DebugMode = 'DEBUG';        % Set this to 'DEBUG' for a fast aproximate output, for complete soundfield contruction replace with ''
Resolution = 100;           % The higher the resolution the better the sampling rate and resultant image
NumberOfPlanewaves = [100 100];   % The more planewaves the more accurate the sound field

% Frequencies_in_zones = [ 500 1500 2000 2500 4000; ...
%                          300 1000 2200 3000 3500];
% Frequency_weights    = [ 0.2  1.5  2.5  5.0 10.0; ...
%                          0.1  1.2  2.8  7.0  8.0];  
                     
Frequencies_in_zones = [3000 ;...%4000; ...  % Zone 1 frequencies
                        3000 ];%3000];     % Zone 2 frequencies
Frequency_weights    = [   0.05   ;...%5 ; ...  % Importance weight for quiet zone of Zone 1
                           0.05   ];%5];      % Importance weight for quiet zone of Zone 2
                     
Angles_in_zones      = [15; ...         % Planewave angle Zone 1
                        15];           % Planewave angle Zone 2
Radius_of_zones      = [0.3; ...        % Zone 1 radius
                        0.3];           % Zone 2 radius
Position_of_zones    = [0.6, 180; ...   % Zone 1 distance and angle from origin
                        0.6, 0];       % Zone 2 distance and angle from origin

Reproduction_Radius = 1.0; % Metres                    

[zones, frequencies] = size(Frequencies_in_zones);  % Currently only support for two zones using Orthogonal Basis Expansion

%zones = 1; % Force just one zone to be created. (Comment out for creation of both zones)

sf          = [];
Soundfield  = [];
BrightField = [];
QuietField  = [];
Bright_d    = [];

%% Calculation & Production of Soundfield

for z = 1:zones
    fprintf('\n====== Wideband Zone %d ======\n\n', z);
    for f = 1:frequencies
        %%
        quiet      = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones(z, f), 0, Radius_of_zones(3-z), 'quiet');
        bright     = Orthogonal_Basis_Expansion.spatial_zone( Frequencies_in_zones(z, f), 0, Radius_of_zones(z),   'pw', 1.0, Angles_in_zones(z));
        quiet.res  = Resolution;
        bright.res = Resolution;
        quiet      = quiet.setDesiredSoundfield(true, 'suppress_output');
        if strcmp(DebugMode, 'DEBUG')
            bright = bright.setDesiredSoundfield();
        else
            bright = bright.setDesiredSoundfield(true, 'suppress_output');
        end
        
        %%
        sf = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
        sf = sf.addSpatialZone(quiet,  Position_of_zones(3-z,1), Position_of_zones(3-z,2));
        sf = sf.addSpatialZone(bright, Position_of_zones(z,1),   Position_of_zones(z,2));
        
        %%
        sf.QuietZ_Weight = Frequency_weights(z, f);
        sf = sf.setN(NumberOfPlanewaves(z));
        sf = sf.createSoundfield(DebugMode, Reproduction_Radius);
        sf = sf.norm_soundfield;        
        
        if (f==1) && (z==1)
            Soundfield = sf.Soundfield_desired;
            BrightField = sf.Bright_Field;
            QuietField = sf.Quiet_Field;
            Bright_d = sf.Bright_Zone.Soundfield_d .* sf.Bright_Zone.Soundfield_d_mask;
        else
            Soundfield = Soundfield + sf.Soundfield_desired;
            if (z==1)
                BrightField = BrightField + sf.Bright_Field;
                QuietField  = QuietField  +  sf.Quiet_Field;
            elseif (z==2)                
                BrightField = BrightField + sf.Quiet_Field;
                QuietField  = QuietField + sf.Bright_Field;
            end
        end

    end
end
Bright_Error = mag2db( sum(sum( abs(real(Bright_d) - real(BrightField)) .^ 2)) / ...
                       sum(sum( abs(         real(Bright_d)           ) .^ 2)) )
Quiet_Error  = mag2db( sum(sum( abs(real(Bright_d) - real(QuietField)) .^ 2)) / ...
                       sum(sum( abs(         real(Bright_d)           ) .^ 2)) )


%% Analysis of Resultant Soundfield
figure(1)
sf.plotSoundfield(Soundfield);
figure(2)
sf.plotSoundfield(BrightField);
figure(3)
sf.plotSoundfield(QuietField);
% figure(2)
% sf.plotSoundfield(BrightField, gray);
% 
% rotated_bright = imrotate(BrightField,Angles_in_zones(1));
% figure(3)
% sf.plotSoundfield(rotated_bright, gray);
% 
% sample = mean(real(rotated_bright));
% sample = sample(sample ~= 0); % TODO: Becareful of this !!!
% figure(4)
% stem(sample);







%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
