clc;
clear;
%close all;

%% Info
Version = {'0weight_withFlatMask', '10000weight_withFlatMask', 'Privacy_Weighted_FlatMask', 'PrivacyWeighted_v5'};
Titles  = {'No Zone Weight with White Noise Masker';
    'Large Uniform Zone Weight with White Noise Masker';
    'Psychoacoustic-Based Zone Weight with White Noise Masker';
    'Psychoacoustic-Based Zone Weight with Zone-Reversed Psychoacoustic-Based Noise Masker'};
LUT_resolution = '512f_256w';
loudspeakers   = 65;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres

figure(2);
for v = 1:length(Version)
    %% Create paths
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = ['+Results\' Output_file_path_ext 'PESQ_Results_' Version{v} '.csv'];
    
    
    % Read results
    [Audio_filename,PESQ_MOS_Bright,PESQ_MOS_Quiet] = Results.import_PESQ(Results_filepath);
    
    NoiseLevel = - cell2mat(textscan([Audio_filename{:}],'%*[^-]%*[^_]_%*[^_]_%d',length(Audio_filename)));
    
    Num_Audio_Files = max(histcounts(NoiseLevel));
    NoiseLevel_Vec = reshape(NoiseLevel,Num_Audio_Files,size(NoiseLevel,1)/Num_Audio_Files);
    NoiseLevel_Vec = double(fliplr(NoiseLevel_Vec(1,:)));
    
    %% Fit a trendline to the results
    [Bright_trend, ~] = Results.createFit(double(NoiseLevel),PESQ_MOS_Bright);
    [Quiet_trend, ~]  = Results.createFit(double(NoiseLevel),PESQ_MOS_Quiet);
    
    %% Calculate areas
    PESQ_Bright_border = flipud([min(reshape(PESQ_MOS_Bright,Num_Audio_Files,size(PESQ_MOS_Bright,1)/Num_Audio_Files)); ...
        max(reshape(PESQ_MOS_Bright,Num_Audio_Files,size(PESQ_MOS_Bright,1)/Num_Audio_Files))]');
    
    PESQ_Bright_area = [PESQ_Bright_border(:,1), PESQ_Bright_border(:,2) - PESQ_Bright_border(:,1)];
    
    PESQ_Quiet_border = flipud([min(reshape(PESQ_MOS_Quiet,Num_Audio_Files,size(PESQ_MOS_Quiet,1)/Num_Audio_Files)); ...
        max(reshape(PESQ_MOS_Quiet,Num_Audio_Files,size(PESQ_MOS_Quiet,1)/Num_Audio_Files))]');
    
    PESQ_Quiet_area = [PESQ_Quiet_border(:,1), PESQ_Quiet_border(:,2) - PESQ_Quiet_border(:,1)];
    
    %% Plot results
    subplot(length(Version)/2,2,v);
    %a_bright = area(NoiseLevel_Vec, PESQ_Bright_area,'LineStyle','none'); hold on;
    
    %a_quiet  = area(NoiseLevel_Vec, PESQ_Quiet_area,'LineStyle','none');
    
    scatter(NoiseLevel,PESQ_MOS_Bright,'bo'); hold on;
    scatter(NoiseLevel,PESQ_MOS_Quiet,'rx');
    
    plot(Bright_trend,'b-');
    plot(Quiet_trend,'r--');
    
    hold off;
    
    title(Titles{v});
    axis([-50 5 0 4.7]);
    grid on; legend({'Bright Zone';'Quiet Zone'});
    
    ylabel('PESQ (MOS)');
    xlabel({'Noise Mask (dB)'; '(with reference to Quiet Zone)'});
    %a_bright(1).FaceColor= 'none';
    %a_bright(2).FaceColor= [0.9 0.9 1];
    
    %a_quiet(1).FaceColor= 'none';
    %a_quiet(2).FaceColor= [1 0.9 0.9];
end