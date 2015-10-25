clc;
%clear;
close all;


%%
ColorInd = 4; %this is for whatever the PESQ colour will be plotted as

%% Info
%LUT_resolution = '512f_256w';
LUT_resolution = '512f_32w';

% loudspeakers   = 295;  % Number of loudspeakers
% speaker_arc    = 360;  % Degrees
%loudspeakers   = 32;  % Number of loudspeakers
speaker_arc    = 180;  % Degrees
loudspeakers   = 16;  % Number of loudspeakers
% speaker_arc    = 180;  % Degrees

speaker_radius = 1.5; % Metres
Num_Receivers = 32;

pw_angle = {0; ...
            0; ...
            15; ...
            15; ...
            90; ...
            90;};
pw_angle = {15};

mask_type_ = {'Flat Mask'; ...
             'Zone Weighted Mask'; ...
             'Flat Mask'; ...
             'Zone Weighted Mask'; ...
             'Flat Mask'; ...
             'Zone Weighted Mask';};
mask_type_ = {'Zone Weighted Mask Alias Ctrl';};

for i=1:length(mask_type_)
    mask_type{1,i} = strrep(strrep(mask_type_{i},'Weighted','Weight'),' ','');
end
         

% % ROOM 1
% % Anechoic6
Room_Size = [10 10 10]; %Anechoic
Wall_Absorption_Coeff = 1.0;

% % ROOM 2
% % Small Open Plan Office
% Room_Size = [4 9 3];   %Small Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 3
% % Medium Open Plan Office
% Room_Size = [8 10 3];   %Medium Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 4
% % Cafe / Restaurant
%Room_Size = [9 14 3];   %Cafe/Restaurant
%Wall_Absorption_Coeff = 0.3;


Reproduction_Centre = Room_Size / 2; %[5 5 5];
if Wall_Absorption_Coeff >= 1.0
    room_type = 'Anechoic';
else
    room_type = 'Reverberant';
end
%{
Version = {['10000weight__with' mask_type{1}]; ...
           ['10000weight__with' mask_type{2}]; ...
           ['10000weight__with' mask_type{3}]; ...
           ['10000weight__with' mask_type{4}]; ...
           ['10000weight__with' mask_type{5}]; ...
           ['10000weight__with' mask_type{6}]};
  %}     
Version = {['10000weight__with' mask_type{1}]};
    %{   
Titles  = {['' mask_type_{1} ]; ... ' and \theta=' num2str(pw_angle{1}) '°']; ...
           ['' mask_type_{2} ]; ... ' and \theta=' num2str(pw_angle{2}) '°']; ...
           ['' mask_type_{3} ]; ... ' and \theta=' num2str(pw_angle{3}) '°']; ...
           ['' mask_type_{4} ]; ... ' and \theta=' num2str(pw_angle{4}) '°']; ...
           ['' mask_type_{5} ]; ... ' and \theta=' num2str(pw_angle{5}) '°']; ...
           ['' mask_type_{6} ]};% ' and \theta=' num2str(pw_angle{6}) '°']};
%}
Titles  = {['' mask_type_{1} ];};% ' and \theta=' num2str(pw_angle{1}) '°']; ...


%Figure Output Settings
DocumentPath = 'tex\latex\Intelligibility';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9)/10;% + 6.35 + 88.9)/10; %IEEE full text width
aspect_ratio = 3/3;
FontSize = 9;
Font = 'Times';
lineWid = 0.1;


%%
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

%%

h = figure(1);
for v = 1:length(Version)
    %% Create paths
    ResultsPath = ['+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = [ResultsPath Output_file_path_ext 'SI_Results_' num2str(pw_angle{v}) 'pwAngle_' Version{v} '.csv'];
    
    
% Read results
    [NoiseLevel,SI_WC_Bright,SI_WC_Quiet,ConfInt_Bright_Low,ConfInt_Bright_Up,ConfInt_Quiet_Low,ConfInt_Quiet_Up] = Results.import_SpeechIntelligibility_Reverb(Results_filepath);
        
    
% Generate Plottable data matrices and vectors
    [Hrz_Vec, Res_Bright_Matrix, Res_Bright_trend, Res_Bright_area, Res_Bright_CI, CI_vec] = ...
        Results.generatePlotData( NoiseLevel, SI_WC_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up, 'smoothingspline', [2 2]);
    [~, Res_Quiet_Matrix, Res_Quiet_trend, Res_Quiet_area, Res_Quiet_CI, ~] = ...
        Results.generatePlotData( NoiseLevel, SI_WC_Quiet, ConfInt_Quiet_Low, ConfInt_Quiet_Up, 'smoothingspline', [2 2]);
    
    
    %% Plot results
    for pl = 1:1 %pl = 1:2
        if pl == 1
            h=figure(1);
%             if mod(length(Version),2)~=1
%                 subplot(length(Version)/2,2,v);
%             else
%                 subplot((length(Version)+1)/2,2,v);
%             end
            if length(Version) ~= 1
                subplot(length(Version)/2,2,v);
            end
        else
            h_sub(v)=figure(100 + v);
        end
        set(gca,'ColorOrderIndex',1);
        arB = area(CI_vec, Res_Bright_area,'LineStyle','none'); hold on;
        arQ  = area(CI_vec, Res_Quiet_area,'LineStyle','none');        
        
        set(gca,'ColorOrderIndex',1);
        erB = errorbar(Hrz_Vec,mean(Res_Bright_Matrix),Res_Bright_CI(:,1),Res_Bright_CI(:,2),...
            'o', 'LineWidth',lineWid, 'MarkerSize',4);
        erQ = errorbar(Hrz_Vec,mean(Res_Quiet_Matrix),Res_Quiet_CI(:,1),Res_Quiet_CI(:,2),...
            'r^', 'LineWidth',lineWid, 'MarkerSize',4);
        
        set(gca,'ColorOrderIndex',1);
        plB = plot(Res_Bright_trend,'-'); set(plB,'LineWidth',lineWid);
        plQ = plot(Res_Quiet_trend,'r-'); set(plQ,'LineWidth',lineWid);
        
        if mod(v,2) == 0
            txt_pos = [-38, 55];
        else
            txt_pos = [0, 55];
        end
       % if v~=1
            text(txt_pos(1),txt_pos(2),['\theta = ' num2str(pw_angle{v}) '°'],'fontname',Font,'fontsize',FontSize)
%         else
%             text(-8,85,['\theta = ' num2str(pw_angle{v}) '°'],'fontname',Font,'fontsize',FontSize)
        %end               
        
        axis([min(Hrz_Vec)-2 20+2 0 105]);
        grid on; grid minor;
        
        legend off;
        
        set(gca,'YTick', round(linspace(0,100,6),1) );
        if mod(v,2)~=0
            ylabel({'STOI (%WC)'});
        else            
            ylabel('');           
        end
        if mod(v,2)==0
             set(gca,'YTickLabel',[]);
        end
        
        if v==(length(Version)-0) || v==(length(Version)-1)
            xlabel({'Noise Mask ({\it{G}}_{dB})'});
        else
            xlabel('');
            set(gca,'XTickLabel',[]);
        end;
        
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
        
        
        
        ax_pos = get(gca,'Position');
        if mod(v,2) == 0
            set(gca,'Position', ax_pos + [-0.08 0 0 0]);
        end
        
        ax_pos = get(gca,'Position');
        if v ~= 1 && v ~= 2
            set(gca,'Position', ax_pos + [0 0.07*round((v-2)/2,0) 0 0]);
        end
        
        set(gca,    'FontSize',FontSize-1, ...
            'FontName',Font);
        
        if v==1 || v==2
            title(Titles{v}, ...
                'FontSize',FontSize, ...
                'FontName',Font);
        end

        
        hold off;
        
    end
    G_=-40:0.1:20;
    contrast_{v}=[Res_Bright_trend(G_)-Res_Quiet_trend(G_)];

end
set(gcf, 'Units','centimeters', 'Color','w');
fig_pos = get(gcf,'Position');
set(gcf, 'PaperUnits','centimeters', ...
    'Position', [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio], ...
    'PaperSize', [plot_width plot_width/aspect_ratio]);

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