clc;
clear;
close all;

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



%Figure Output Settings
DocumentPath = 'tex\latex\Intelligibility';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10);% + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 9/3;
FontSize = 9;
Font = 'Times';


%%

h = figure(1);
for v = 1:length(Version)
    %% Create paths
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = ['Z:\+Results\' Output_file_path_ext 'SI_Results_' Version{v} '.csv'];
    
    
    % Read results
    [Audio_filename,SI_WC_Bright,SI_WC_Quiet] = Results.import_SpeechIntelligibility(Results_filepath);
    
    NoiseLevel = - cell2mat(textscan([Audio_filename{:}],'%*[^-]%*[^_]_%*[^_]_%d',length(Audio_filename)));
    
    Num_Audio_Files = max(histcounts(NoiseLevel));
    
    SI_WC_Bright_Matrix = fliplr(reshape(SI_WC_Bright,Num_Audio_Files,size(SI_WC_Bright,1)/Num_Audio_Files));
    SI_WC_Quiet_Matrix  = fliplr(reshape(SI_WC_Quiet ,Num_Audio_Files,size(SI_WC_Quiet ,1)/Num_Audio_Files));
    NoiseLevel_Vec = reshape(NoiseLevel,Num_Audio_Files,size(NoiseLevel,1)/Num_Audio_Files);
    NoiseLevel_Vec = double(fliplr(NoiseLevel_Vec(1,:)));
    
        %% Calculate confidence intervals
    SI_WC_Bright_CI = Tools.confidence_intervals(SI_WC_Bright_Matrix, 95);
    SI_WC_Quiet_CI = Tools.confidence_intervals(SI_WC_Quiet_Matrix, 95);
    
    %% Fit a trendline to the results
    Fit_Options = fitoptions( 'Method', 'SmoothingSpline' );
    Fit_Options.SmoothingParam = 1.0;
    [Bright_trend, ~] = Results.createFit(double(NoiseLevel),SI_WC_Bright, 'smoothingspline', Fit_Options);
    [Quiet_trend, ~]  = Results.createFit(double(NoiseLevel),SI_WC_Quiet, 'smoothingspline', Fit_Options);
    
    temp1 = repmat(SI_WC_Bright_CI(:,1)',Num_Audio_Files,1); temp1 = temp1(:);
    temp2 = repmat(SI_WC_Bright_CI(:,2)',Num_Audio_Files,1); temp2 = temp2(:);
    Bright_CI_trend = struct('Upper', ...
                     Results.createFit(double(NoiseLevel),temp1,'smoothingspline',Fit_Options), ...
                     'Lower', ...
                     Results.createFit(double(NoiseLevel),temp2,'smoothingspline',Fit_Options));
                 
    temp1 = repmat(SI_WC_Quiet_CI(:,1)',Num_Audio_Files,1); temp1 = temp1(:);
    temp2 = repmat(SI_WC_Quiet_CI(:,2)',Num_Audio_Files,1); temp2 = temp2(:);
    Quiet_CI_trend = struct('Upper', ...
                     Results.createFit(double(NoiseLevel),temp1,'smoothingspline',Fit_Options), ...
                     'Lower', ...
                     Results.createFit(double(NoiseLevel),temp2,'smoothingspline',Fit_Options));
                 
    %Trendline plot vectors
    CI_x= linspace(NoiseLevel_Vec(1),NoiseLevel_Vec(end),100);
    
    %% Calculate areas
    %SI_WC_Bright_area = [mean(SI_WC_Bright_Matrix)' + SI_WC_Bright_CI(:,1), SI_WC_Bright_CI(:,2) - SI_WC_Bright_CI(:,1)];
    SI_WC_Bright_area = flip([flip(Bright_trend(CI_x)) + Bright_CI_trend.Upper(CI_x), ...
                         Bright_CI_trend.Lower(CI_x) - Bright_CI_trend.Upper(CI_x)]);
    
    %SI_WC_Quiet_area = [mean(SI_WC_Quiet_Matrix)' + SI_WC_Quiet_CI(:,1), SI_WC_Quiet_CI(:,2) - SI_WC_Quiet_CI(:,1)];
    SI_WC_Quiet_area = flip([flip(Quiet_trend(CI_x)) + Quiet_CI_trend.Upper(CI_x), ...
                        Quiet_CI_trend.Lower(CI_x) - Quiet_CI_trend.Upper(CI_x)]);
    
    
    %% Plot results
    for pl = 1:2
        if pl == 1
            h=figure(1);
            subplot(length(Version)/2,2,v);
        else
            h_sub(v)=figure(100 + v);
        end
        arB = area(CI_x, SI_WC_Bright_area,'LineStyle','none'); hold on;
        arQ  = area(CI_x, SI_WC_Quiet_area,'LineStyle','none');
        
        %scB = scatter(NoiseLevel,SI_WC_Bright,'bo');
        %scQ = scatter(NoiseLevel,SI_WC_Quiet,'rx');
        set(gca,'ColorOrderIndex',1);
        erB = errorbar(NoiseLevel_Vec,mean(SI_WC_Bright_Matrix),SI_WC_Bright_CI(:,1),SI_WC_Bright_CI(:,2),...
            'o');
        erQ = errorbar(NoiseLevel_Vec,mean(SI_WC_Quiet_Matrix),SI_WC_Quiet_CI(:,1),SI_WC_Quiet_CI(:,2),...
            'r^');
        
        set(gca,'ColorOrderIndex',1);
        plB = plot(Bright_trend,'-');
        plQ = plot(Quiet_trend,'r--');
        
        hold off;
        
        title(Titles{v});
        axis([-46 1 0 100]);
        grid on;
        legend([erB, erQ], {'Bright Zone';'Quiet Zone'},'Location','southwest');
        
        ylabel({'Speech Intelligibility';'(% Words Correct)'});
        xlabel({'Noise Mask (dB)'});%; '(with reference to Quiet Zone)'});
        
        arB(1).FaceColor= 'none';
        arB(2).FaceColor= [0.8 0.9 1];
        
        arQ(1).FaceColor= 'none';
        arQ(2).FaceColor= [1 0.9 0.9];
    end
end

%% 
figure(100 + 2);
set(gca,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1]);
grid on;box on; grid minor;
        
%t(1).Position(2) = t(1).Position(2);
factor1 = 2.3;
factor2 = 0.05;
offset_cm = FontSize*254/720 /10 * factor1;
set(gcf, 'PaperUnits','centimeters', ...
    'PaperSize', [plot_width, plot_width/aspect_ratio + offset_cm], ...
    'PaperPosition',[0 factor2 plot_width,plot_width/aspect_ratio + offset_cm]);

t(1) = title('Contrast in Speech Intelligibility from Added Noise', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
t(1).Position(2) = t(1).Position(2) + FontSize;

le = legend();
%set(le,'Box','off');
le.Position = le.Position + [0.10 -0.020 0 0];

%print(['-d' print_fmt], [DocumentPath '\SI_Contrast_Enhancement_All']);

% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;

%close all;