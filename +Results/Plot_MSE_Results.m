clc;
clear;
close all;

%% Info
Weights_Vec = [1e-2 1e0 1e2 1e4];
Version = {'0weight', '1weight', '100weight', '10000weight'};
Titles  = Version;
LUT_resolution = '512f_256w';
loudspeakers   = 65;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres



%Figure Output Settings
DocumentPath = 'tex\latex\APSIPA';
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
    Results_filepath = ['Z:\+Results\' Output_file_path_ext 'MSE_Results_' Version{v} '.csv'];
    
    
    % Read results
    [Audio_filename,MSE_Bright_Matrix(:,v),MSE_Quiet_Matrix(:,v)] = Results.import_MSE(Results_filepath); %#ok<SAGROW>
    
end
    Num_Audio_Files = length(Audio_filename);

    % Convert magnitude to dB
    MSE_Bright_Matrix = mag2db(MSE_Bright_Matrix);
    %MSE_Quiet_Matrix = mag2db(MSE_Quiet_Matrix);
    
    MSE_Bright = MSE_Bright_Matrix(:);
    %MSE_Quiet = MSE_Quiet_Matrix(:);
    Weights = repmat(Weights_Vec,Num_Audio_Files,1);
    Weights = Weights(:);

    
    %% Calculate confidence intervals
    MSE_Bright_CI = Tools.confidence_intervals(MSE_Bright_Matrix, 95);
    %MSE_Quiet_CI = Tools.confidence_intervals(MSE_Quiet_Matrix, 95);
    
    %% Fit a trendline to the results    
    Fit_Options = fitoptions( 'Method', 'SmoothingSpline' );
    Fit_Options.SmoothingParam = 1.0;
    [Bright_trend, ~] = Results.createFit(double(log10(Weights)+2),MSE_Bright,'smoothingspline',Fit_Options);
    %[Quiet_trend, ~]  = Results.createFit(double(Weights),MSE_Quiet);
    
    temp1 = repmat(MSE_Bright_CI(:,1)',Num_Audio_Files,1); temp1 = temp1(:);
    temp2 = repmat(MSE_Bright_CI(:,2)',Num_Audio_Files,1); temp2 = temp2(:);
    Bright_CI_trend = struct('Upper', ...
                     Results.createFit(double(log10(Weights)+2),temp1,'smoothingspline',Fit_Options), ...
                     'Lower', ...
                     Results.createFit(double(log10(Weights)+2),temp2,'smoothingspline',Fit_Options));
    %Trendline plot vectors
    x = logspace(log10(1e-2),log10(1e4),100);
    xx= linspace(-2,4,100)+2;
    
    %% Calculate areas
    MSE_Bright_area = [Bright_trend(xx) + Bright_CI_trend.Upper(xx), Bright_CI_trend.Lower(xx) - Bright_CI_trend.Upper(xx)];
    %MSE_Quiet_area = [mean(MSE_Quiet_Matrix)' + MSE_Quiet_CI(:,1), MSE_Quiet_CI(:,2) - MSE_Quiet_CI(:,1)];
    
    
    %% Plot results
    
    arB = area(x, MSE_Bright_area,'LineStyle','none'); hold on;
    %arQ  = area(Weights_Vec, MSE_Quiet_area,'LineStyle','none');
    
    set(gca,'ColorOrderIndex',1);
    erB = errorbar(Weights_Vec,mean(MSE_Bright_Matrix),MSE_Bright_CI(:,1),MSE_Bright_CI(:,2),'o');
    %erQ = errorbar(Weights_Vec,mean(MSE_Quiet_Matrix),MSE_Quiet_CI(:,1),MSE_Quiet_CI(:,2),'rx');
    
    set(gca,'ColorOrderIndex',1);
    plB = plot(x, Bright_trend(xx),'-');
    %plQ = plot(Quiet_trend,'r--');
    
    hold off;
    
    title(Titles{v});
    xlim([1*10^(-2.5) 1*10^(4.5)]);
    ylim([-85 -65]);
    set(gca, 'XTick', [1e-2, 1e0, 1e2, 1e4]);
    set(gca,'XScale','log');
    %legend([arB(2), arQ(2)], {'Bright Zone';'Quiet Zone'});
    %legend([scB, scQ], {'Bright Zone';'Quiet Zone'});
    legend([erB], ... %, plQ], 
        {'Bright Zone'}, ...%;'Quiet Zone'},
        'Location','southeast');
    
    ylabel('MSE (dB)');
    xlabel('Weight');
    
    arB(1).FaceColor= 'none';
    arB(2).FaceColor= [0.8 0.9 1];
    
    arQ(1).FaceColor= 'none';
    arQ(2).FaceColor= [1 0.9 0.9];



%% 
grid on; box on;
set(gca,    'FontSize',FontSize, ...
            'FontName',Font, ...
            'PlotBoxAspectRatio',[1,1/aspect_ratio,1]);

        
%t(1).Position(2) = t(1).Position(2);
factor1 = 2.5;
factor2 = factor1*2;
offset_cm = FontSize*254/720 /10 * factor1;
set(gcf, 'PaperUnits','centimeters', ...
    'PaperSize', [plot_width, plot_width/aspect_ratio + offset_cm], ...
    'PaperPosition',[0 offset_cm/factor2 plot_width,plot_width/aspect_ratio + offset_cm]);

t(1) = title('MSE of Weighted Reproduced Speech', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
%t(1).Position(2) = t(1).Position(2) + FontSize/2;

print(['-d' print_fmt], [DocumentPath '\Weighting_Affect_on_MSE']);

% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;

