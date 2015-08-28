clc;
clear;
%close all;

%% Add to existing figure
add2fig = 1;
ColorOrderIndex = 4;

%% Info
result_type = 'PESQ';
LUT_resolution = '512f_256w';
loudspeakers   = 295;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres
Num_Receivers = 32;

pw_angle = {0; ...
            15; ...
            15; ...
            90; ...
            90;};
mask_type = {'FlatMask'; ...
             'FlatMask'; ...
             'ZoneWeightMask'; ...
             'FlatMask'; ...
             'ZoneWeightMask';};
         
Room_Size = [10 10 10];% ...
             %[5 6 ]; ...
             %[5 6 4]; ...
             %[5 6 4];};
Reproduction_Centre = [5 5 5];% ...
                      % [2.5 3 ]; ...
                      % [2.5 3 2]; ...
                      % [2.5 3 2];};
Wall_Absorption_Coeff = 1.0;%{1.0; 0.3; 1.0; 0.3;};% 0.3};


Version = {['10000weight__with' mask_type{1}]; ...
           ['10000weight__with' mask_type{2}]; ...
           ['10000weight__with' mask_type{3}]; ...
           ['10000weight__with' mask_type{4}]; ...
           ['10000weight__with' mask_type{5}]};
       
Titles  = {['Large Uniform Zone Weight with ' mask_type{1} ' - Anechoic ' num2str(pw_angle{1}) 'deg']; ...
           ['Large Uniform Zone Weight with ' mask_type{2} ' - Anechoic ' num2str(pw_angle{2}) 'deg']; ...
           ['Large Uniform Zone Weight with ' mask_type{3} ' - Anechoic ' num2str(pw_angle{3}) 'deg']; ...
           ['Large Uniform Zone Weight with ' mask_type{4} ' - Anechoic ' num2str(pw_angle{4}) 'deg']; ...
           ['Large Uniform Zone Weight with ' mask_type{5} ' - Anechoic ' num2str(pw_angle{5}) 'deg']};


%Figure Output Settings
DocumentPath = 'tex\latex\Intelligibility';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10);% + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 9/3;
FontSize = 9;
Font = 'Times';


%%
room = strrep(sprintf(strrep(repmat('%d',1,length(Room_Size)),'d%','d %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%d',1,length(Reproduction_Centre)),'d%','d %'),Reproduction_Centre),' ','x');

%%

h = figure(1);
for v = 1:length(Version)
    %% Create paths
    ResultsPath = ['+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = [ResultsPath Output_file_path_ext result_type '_Results_' num2str(pw_angle{v}) 'pwAngle_' Version{v} '.csv'];
    
    
% Read results
    [NoiseLevel,Result_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] = Results.import_PESQ_Reverb(Results_filepath);
        
    Num_Audio_Files = max(histcounts(NoiseLevel));
    
% Create Matrices of Results
    %Masking noise Level
    NoiseLevel_Matrix = fliplr(reshape(NoiseLevel,size(NoiseLevel,1)/Num_Audio_Files,Num_Audio_Files)');
    NoiseLevel_Vec = double(NoiseLevel_Matrix(1,:));
    %Speech Intelligibility Results
    Result_Bright_Matrix = fliplr(reshape(Result_Bright,size(Result_Bright,1)/Num_Audio_Files,Num_Audio_Files)');
    % Confidence Intervals from Spatial Sampling Points
    ConfInt_Bright_Low_M = fliplr(reshape(ConfInt_Bright_Low,size(ConfInt_Bright_Low,1)/Num_Audio_Files,Num_Audio_Files)');
    ConfInt_Bright_Up_M  = fliplr(reshape(ConfInt_Bright_Up ,size(ConfInt_Bright_Up ,1)/Num_Audio_Files,Num_Audio_Files)');
    
% Order Matrices of Results
    [NoiseLevel_Vec,I] = sort(NoiseLevel_Vec);
    %Masking noise Level
    NoiseLevel_Matrix   = NoiseLevel_Matrix(:,I);
    NoiseLevel = flip(NoiseLevel_Matrix(:));
    %Speech Intelligibility Results
    Result_Bright_Matrix = Result_Bright_Matrix(:,I);
    % Confidence Intervals from Spatial Sampling Points
    ConfInt_Bright_Low_M = ConfInt_Bright_Low_M(:,I);
    ConfInt_Bright_Up_M  = ConfInt_Bright_Up_M(:,I);
    
    %Average Spatial Sampling Confidence Intervals
    ConfInt_Bright_M = [mean(ConfInt_Bright_Low_M,1)' mean(ConfInt_Bright_Up_M,1)'];
    
        %% Calculate confidence intervals
    Result_Bright_CI = Tools.confidence_intervals(Result_Bright_Matrix, 95);
    
    %% Fit a trendline to the results
    Fit_Options = fitoptions( 'Method', 'SmoothingSpline' );
    Fit_Options.SmoothingParam = 1.0;
    [Bright_trend, ~] = Results.createFit(double(NoiseLevel),flip(Result_Bright_Matrix(:)), 'smoothingspline', Fit_Options);
    
    temp1 = repmat(Result_Bright_CI(:,1)',Num_Audio_Files,1); temp1 = temp1(:);
    temp2 = repmat(Result_Bright_CI(:,2)',Num_Audio_Files,1); temp2 = temp2(:);
    Bright_CI_trend = struct('Upper', ...
                     Results.createFit(double(NoiseLevel),temp1,'smoothingspline',Fit_Options), ...
                     'Lower', ...
                     Results.createFit(double(NoiseLevel),temp2,'smoothingspline',Fit_Options));
                 
    %Trendline plot vectors
    CI_x= linspace(NoiseLevel_Vec(1),NoiseLevel_Vec(end),100);
    
    %% Calculate areas
    %SI_WC_Bright_area = [mean(SI_WC_Bright_Matrix)' + SI_WC_Bright_CI(:,1), SI_WC_Bright_CI(:,2) - SI_WC_Bright_CI(:,1)];
    Result_Bright_area = flip([flip(Bright_trend(CI_x)) + Bright_CI_trend.Upper(CI_x), ...
                         Bright_CI_trend.Lower(CI_x) - Bright_CI_trend.Upper(CI_x)]);  
    
    %% Plot results
    for pl = 1:2
        if pl == 1
            if add2fig
                h=figure(add2fig);
                ax = h.Children(end - (v-1)*2);
                gca = axes('Position',ax.Position,...
                    'YAxisLocation','right',...
                    'Color','none');
                hold on;
            else
                h=figure(2);
                %if mod(length(Version),2)~=1
                %	subplot(length(Version)/2,2,v);
                %else
                %	subplot((length(Version)+1)/2,2,v);
                %end
                subplot(length(Version),1,v);
            end
        else
            h_sub(v)=figure(100 + v);
        end
        arB = area(CI_x, Result_Bright_area,'LineStyle','none'); hold on;
        
        %scB = scatter(NoiseLevel,SI_WC_Bright,'bo');
        set(gca,'ColorOrderIndex',ColorOrderIndex);
        erB = errorbar(NoiseLevel_Vec,mean(Result_Bright_Matrix),Result_Bright_CI(:,1),Result_Bright_CI(:,2),...
            's');
        
        set(gca,'ColorOrderIndex',ColorOrderIndex);
        plB = plot(Bright_trend,'-');
        
        hold off;
        
        title(Titles{v});
        axis([-41 41 1 4.6]);
        grid on;
        legend(erB, {'Bright Zone'},'Location','northeast');
        
        ylabel({'PESQ (MOS)'});
        xlabel({'Noise Mask (dB)'});%; '(with reference to Quiet Zone)'});
        
        arB(1).FaceColor = 'none';
        arB(2).FaceColor = plB.Color; 
        drawnow; pause(0.05);  % this is important for transparency!
        arB(2).Face.ColorType = 'truecoloralpha';
        arB(2).Face.ColorData(4) = 0.2*255;              
        
        if add2fig
            ylim( [1 4.75*(104/100)] );
            gca.YTick=linspace(1,4.75,6);
            gca.XTick=[];
            gca.XLabel=[];
            gca.Title=[];
        end
        
    end
end

%% 
% figure(100 + 2);
% set(gca,    'FontSize',FontSize, ...
%             'FontName',Font, ...
%             'PlotBoxAspectRatio',[1,1/aspect_ratio,1]);
% grid on;box on; grid minor;
%         
% %t(1).Position(2) = t(1).Position(2);
% factor1 = 2.3;
% factor2 = 0.05;
% offset_cm = FontSize*254/720 /10 * factor1;
% set(gcf, 'PaperUnits','centimeters', ...
%     'PaperSize', [plot_width, plot_width/aspect_ratio + offset_cm], ...
%     'PaperPosition',[0 factor2 plot_width,plot_width/aspect_ratio + offset_cm]);
% 
% t(1) = title('Contrast in Speech Intelligibility from Added Noise', 'FontWeight','bold', 'FontSize',FontSize+1, 'FontName',Font);
% t(1).Position(2) = t(1).Position(2) + FontSize;
% 
% le = legend();
% %set(le,'Box','off');
% le.Position = le.Position + [0.10 -0.020 0 0];
% 
% print(['-d' print_fmt], [DocumentPath '\SI_Contrast_Enhancement']);
% 
% close all;
% % Update latex File Name DataBase
% Tools.MiKTeX_FNDB_Refresh;
% 
% %close all;