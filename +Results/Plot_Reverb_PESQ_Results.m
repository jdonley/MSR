clc;
clear;
%close all;

%% Add to existing figure
add2fig = 0;
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
         

% % ROOM 1
% % Anechoic6
Room_Size = [10 10 10]; %Anechoic
Wall_Absorption_Coeff = 1.0;

% % ROOM 2
% % Small Open Plan Office
% Room_Size = [4 9 3];   %Small Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 3
% Medium Open Plan Office
% Room_Size = [8 10 3];   %Medium Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 4
% % Cafe / Restaurant
%Room_Size = [9 14 3];   %Cafe/Restaurant
%Wall_Absorption_Coeff = 0.3;


Reproduction_Centre = Room_Size / 2; %[5 5 5];
if Wall_Absorption_Coeff == 1.0
    room_type = 'Anechoic';
elseif Wall_Absorption_Coeff == 0.3
    room_type = 'Reverberant';
end    

Version = {['10000weight__with' mask_type{1}]; ...
           ['10000weight__with' mask_type{2}]; ...
           ['10000weight__with' mask_type{3}]; ...
           ['10000weight__with' mask_type{4}]; ...
           ['10000weight__with' mask_type{5}]};
       
Titles  = {['Large Zone Weight with ' mask_type{1} ', ' room_type ', ' num2str(pw_angle{1}) 'deg']; ...
           ['Large Zone Weight with ' mask_type{2} ', ' room_type ', ' num2str(pw_angle{2}) 'deg']; ...
           ['Large Zone Weight with ' mask_type{3} ', ' room_type ', ' num2str(pw_angle{3}) 'deg']; ...
           ['Large Zone Weight with ' mask_type{4} ', ' room_type ', ' num2str(pw_angle{4}) 'deg']; ...
           ['Large Zone Weight with ' mask_type{5} ', ' room_type ', ' num2str(pw_angle{5}) 'deg']};


%Figure Output Settings
DocumentPath = 'tex\latex\Intelligibility';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9/10);% + 6.35/10 + 88.9/10); %IEEE full page width
aspect_ratio = 9/3;
FontSize = 9;
Font = 'Times';


%%
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

%%

h = figure(1);
for v = 1:length(Version)
    %% Create paths
    ResultsPath = ['+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = [ResultsPath Output_file_path_ext result_type '_Results_' num2str(pw_angle{v}) 'pwAngle_' Version{v} '.csv'];
    
    
% Read results
    [NoiseLevel,Result_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] = Results.import_PESQ_Reverb(Results_filepath);
        
% Generate Plottable data matrices and vectors
    [Hrz_Vec, Res_Bright_Matrix, Res_Bright_trend, Res_Bright_area, Res_Bright_CI, CI_vec] = Results.generatePlotData( NoiseLevel, Result_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up);
    
    
    %% Plot results
    for pl = 1:2
        if pl == 1
            if add2fig
                h=figure(add2fig);
                if ~numel(h.Children)
                    error('No figure exists to add plots too.\nPlease change the ''add2fig'' value to zero.',[]);
                end
                ax = h.Children(end - (v-1)*2);
                ax2 = axes('Position',ax.Position,...
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
        arB = area(CI_vec, Res_Bright_area,'LineStyle','none'); hold on;
        
        %scB = scatter(NoiseLevel,SI_WC_Bright,'bo');
        set(gca,'ColorOrderIndex',ColorOrderIndex);
        erB = errorbar(Hrz_Vec,mean(Res_Bright_Matrix),Res_Bright_CI(:,1),Res_Bright_CI(:,2),...
            's');
        
        set(gca,'ColorOrderIndex',ColorOrderIndex);
        plB = plot(Res_Bright_trend,'-');
        
        hold off;
        
        title(Titles{v});
        axis([min(Hrz_Vec)-1 max(Hrz_Vec)+1 1 4.6]);
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
            %ylim( [1 4.75*(104/100)] );
            %gca.YTick=linspace(1,4.75,6);
            ax2.XTick=[];
            ax2.XLabel=[];
            ax2.Title=[];
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