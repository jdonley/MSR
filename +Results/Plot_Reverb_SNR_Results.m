clc;
clear;
close all;

%% Info
result_type = 'SNR';
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
    [NoiseLevel,Result_Bright,Result_Quiet,ConfInt_Bright_Low,ConfInt_Bright_Up,ConfInt_Quiet_Low,ConfInt_Quiet_Up] = Results.import_SegSNR_Reverb(Results_filepath);
            
% Generate Plottable data matrices and vectors
    [Hrz_Vec, Res_Bright_Matrix, Res_Bright_trend, Res_Bright_area, Res_Bright_CI, CI_vec] = Results.generatePlotData( NoiseLevel, Result_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up);
    [~, Res_Quiet_Matrix, Res_Quiet_trend, Res_Quiet_area, Res_Quiet_CI, ~] = Results.generatePlotData( NoiseLevel, Result_Quiet, ConfInt_Quiet_Low, ConfInt_Quiet_Up);
    
    
    %% Plot results
    for pl = 1:2
        if pl == 1
            h=figure(1);
%             if mod(length(Version),2)~=1
%                 subplot(length(Version)/2,2,v);
%             else
%                 subplot((length(Version)+1)/2,2,v);
%             end
            subplot(length(Version),1,v);
        else
            h_sub(v)=figure(100 + v);
        end
        arB = area(CI_vec, Res_Bright_area,'LineStyle','none'); hold on;
        arQ  = area(CI_vec, Res_Quiet_area,'LineStyle','none');
        
        %scB = scatter(NoiseLevel,SI_WC_Bright,'bo');
        %scQ = scatter(NoiseLevel,SI_WC_Quiet,'rx');
        set(gca,'ColorOrderIndex',1);
        erB = errorbar(Hrz_Vec,mean(Res_Bright_Matrix),Res_Bright_CI(:,1),Res_Bright_CI(:,2),...
            'o');
        erQ = errorbar(Hrz_Vec,mean(Res_Quiet_Matrix),Res_Quiet_CI(:,1),Res_Quiet_CI(:,2),...
            'r^');
        
        set(gca,'ColorOrderIndex',1);
        plB = plot(Res_Bright_trend,'-');
        plQ = plot(Res_Quiet_trend,'r--');
        
        hold off;
        
        title(Titles{v});
        axis([min(Hrz_Vec)-1 max(Hrz_Vec)+1 -20 10]);
        grid on;
        legend([erB, erQ], {'Bright Zone';'Quiet Zone'},'Location','northeast');
        
        ylabel({'Segmental Signal to Noise Ratio (dB)'});
        xlabel({'Noise Mask (dB)'});%; '(with reference to Quiet Zone)'});
        
        
        arB(1).FaceColor= 'none';
        arB(2).FaceColor= plB.Color; 
        drawnow; pause(0.05);  % this is important for transparency!
        arB(2).Face.ColorType = 'truecoloralpha';
        arB(2).Face.ColorData(4) = 0.2*255;
        
        arQ(1).FaceColor= 'none';
        arQ(2).FaceColor= plQ.Color;
        drawnow; pause(0.05);  % this is important for transparency!
        arQ(2).Face.ColorType = 'truecoloralpha';
        arQ(2).Face.ColorData(4) = 0.2*255;
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