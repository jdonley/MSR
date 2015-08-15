clc;
clear;
%close all;

%% Info
Version = {'0weight_withFlatMask', '10000weight_withFlatMask', 'Privacy_Weighted_v6_FlatMask', 'PrivacyWeighted_v5'};
Titles  = {'No Zone Weight with White Noise Masker';
    'Large Uniform Zone Weight with White Noise Masker';
    'Psychoacoustic-Based Zone Weight with White Noise Masker';
    'Psychoacoustic-Based Zone Weight with Zone-Reversed Psychoacoustic-Based Noise Masker'};
LUT_resolution = '512f_256w';
loudspeakers   = {65, 65, 295, 65};  % Number of loudspeakers
%loudspeakers   = {65, 65, 65, 65};  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres

figure(1);
for v = 1:length(Version)
    %% Create paths
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers{v}) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = ['+Results\' Output_file_path_ext 'SI_Results_' Version{v} '.csv'];
    
    
    % Read results
    [Audio_filename,SI_WC_Bright,SI_WC_Quiet] = Results.import_SpeechIntelligibility(Results_filepath);
    
    NoiseLevel = - cell2mat(textscan([Audio_filename{:}],'%*[^-]%*[^_]_%*[^_]_%d',length(Audio_filename)));
    
    Num_Audio_Files = max(histcounts(NoiseLevel));
    
    SI_WC_Bright_Matrix = fliplr(reshape(SI_WC_Bright,Num_Audio_Files,size(SI_WC_Bright,1)/Num_Audio_Files));
    SI_WC_Quiet_Matrix  = fliplr(reshape(SI_WC_Quiet ,Num_Audio_Files,size(SI_WC_Quiet ,1)/Num_Audio_Files));
    NoiseLevel_Vec = reshape(NoiseLevel,Num_Audio_Files,size(NoiseLevel,1)/Num_Audio_Files);
    NoiseLevel_Vec = double(fliplr(NoiseLevel_Vec(1,:)));
    
    %% Fit a trendline to the results
    [Bright_trend, ~] = Results.createFit(double(NoiseLevel),SI_WC_Bright);
    [Quiet_trend, ~]  = Results.createFit(double(NoiseLevel),SI_WC_Quiet);
    NoiseLevel = NoiseLevel_Vec(1):NoiseLevel_Vec(end);
    Contrast_fromTrend = Bright_trend(NoiseLevel) - Quiet_trend(NoiseLevel);
    
    %% Plot results
    figure(1);
    pl = plot(NoiseLevel, Contrast_fromTrend,'LineWidth',2);
    hold on;
    
end

hold off;
title('Speech Intelligibility Contrast for Different Weighting Schemes');
axis([-46 1 0 100]);
grid on;
%legend([arB(2), arQ(2)], {'Bright Zone';'Quiet Zone'});
%legend([scB, scQ], {'Bright Zone';'Quiet Zone'});
legend(Titles,'Location','northwest');

ylabel({'Speech Intelligibility Contrast';'(% Words Correct)'});
xlabel({'Noise Mask (dB)'; '(with reference to Quiet Zone)'});