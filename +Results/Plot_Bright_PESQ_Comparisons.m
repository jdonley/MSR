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

figure(1);
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
    NoiseLevel = NoiseLevel_Vec(1):NoiseLevel_Vec(end);
    
    %% Plot results
    figure(2);
    pl = plot(NoiseLevel, Bright_trend(NoiseLevel),'LineWidth',2);
    hold on;
    
end

hold off;
title('PESQ for Different Weighting Schemes');
axis([-46 1 0 4.6]);
grid on;
%legend([arB(2), arQ(2)], {'Bright Zone';'Quiet Zone'});
%legend([scB, scQ], {'Bright Zone';'Quiet Zone'});
legend(Titles,'Location','northwest');

ylabel('PESQ (MOS)');
xlabel({'Noise Mask (dB)'; '(with reference to Quiet Zone)'});