clc;
clear;
%close all;

%% Add to existing figure
%Plot_ = 'STOI&PESQ';
Plot_ = 'STI&PESQ';
if strcmp(Plot_,'STOI&PESQ')
    Results.Plot_Reverb_STOI_Results;
    Plot_ = 'STOI&PESQ';
elseif strcmp(Plot_,'STI&PESQ')
    Results.Plot_Reverb_STI_Results;
    Plot_ = 'STI&PESQ';
end
add2fig = 1; % If not set to zero but instead set to a postiive number then this script will try and add the plots to an existing figure with that number
ColorInd = 4; % Colour of lines for these plots

%% Info
result_type = 'PESQ';
LUT_resolution = '512f_256w';
loudspeakers   = 295;  % Number of loudspeakers
speaker_arc    = 360;  % Degrees
speaker_radius = 1.5; % Metres
Num_Receivers = 32;

pw_angle = {0; ...
            0; ...
            15; ...
            15; ...
            90; ...
            90;};
mask_type_ = {'Flat Mask'; ...
             'Zone Weighted Mask'; ...
             'Flat Mask'; ...
             'Zone Weighted Mask'; ...
             'Flat Mask'; ...
             'Zone Weighted Mask';};
for i=1:length(mask_type_)
    mask_type{1,i} = strrep(strrep(mask_type_{i},'Weighted','Weight'),' ','');
end

if strcmp(Plot_,'STOI&PESQ')
    rs =[1];
else
    rs = [2 3 4];
end

for r_ = 1:length(rs)
    r = rs(r_);
    if r == 1
        % % ROOM 1
        % % Anechoic6
        Room_Size = [10 10 10]; %Anechoic
        Wall_Absorption_Coeff = 1.0;
    elseif r == 2
        % % ROOM 2
        % % Small Open Plan Office
        Room_Size = [4 9 3];   %Small Open Plan Office
        Wall_Absorption_Coeff = 0.3;
    elseif r == 3
        % % ROOM 3
        % % Medium Open Plan Office
        Room_Size = [8 10 3];   %Medium Open Plan Office
        Wall_Absorption_Coeff = 0.3;
    elseif r == 4
        % % ROOM 4
        % % Cafe / Restaurant
        Room_Size = [9 14 3];   %Cafe/Restaurant
        Wall_Absorption_Coeff = 0.3;
    end


Reproduction_Centre = Room_Size / 2; %[5 5 5];
if Wall_Absorption_Coeff >= 1.0
    room_type = 'Anechoic';
else
    room_type = 'Reverberant';
end

Version = {['10000weight__with' mask_type{1}]; ...
           ['10000weight__with' mask_type{2}]; ...
           ['10000weight__with' mask_type{3}]; ...
           ['10000weight__with' mask_type{4}]; ...
           ['10000weight__with' mask_type{5}]; ...
           ['10000weight__with' mask_type{6}]};
       
Titles  = {['Large Zone Weight with ' mask_type_{1} ', ' room_type ', ' num2str(pw_angle{1}) 'deg']; ...
           ['Large Zone Weight with ' mask_type_{2} ', ' room_type ', ' num2str(pw_angle{2}) 'deg']; ...
           ['Large Zone Weight with ' mask_type_{3} ', ' room_type ', ' num2str(pw_angle{3}) 'deg']; ...
           ['Large Zone Weight with ' mask_type_{4} ', ' room_type ', ' num2str(pw_angle{4}) 'deg']; ...
           ['Large Zone Weight with ' mask_type_{5} ', ' room_type ', ' num2str(pw_angle{5}) 'deg']; ...
           ['Large Zone Weight with ' mask_type_{6} ', ' room_type ', ' num2str(pw_angle{6}) 'deg']};


%Figure Output Settings
DocumentPath = 'tex\latex\ICASSP';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = (88.9)/10;% + 6.35 + 88.9)/10; %IEEE full text width
aspect_ratio = 3/3;
FontSize = 9;
Font = 'Times';
lineWid = 0.5;


%%
room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');

%%

h = figure(1);
drawnow; pause(0.05);  % this is important for transparency!
for v = 1:length(Version)
    %% Create paths
    ResultsPath = ['+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];
    Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
    Results_filepath = [ResultsPath Output_file_path_ext result_type '_Results_' num2str(pw_angle{v}) 'pwAngle_' Version{v} '.csv'];
    
    
% Read results
    [NoiseLevel,Result_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] = Results.import_PESQ_Reverb(Results_filepath);
        
% Generate Plottable data matrices and vectors
    [Hrz_Vec, Res_Bright_Matrix, Res_Bright_trend, Res_Bright_area, Res_Bright_CI, CI_vec] = ...
        Results.generatePlotData( NoiseLevel, Result_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up, 'smoothingspline', [1.8 1.8]);
    
    
    %% Plot results
    for pl = 1:1 %pl = 1:2
        if pl == 1
            if add2fig
                h=figure(add2fig);
                if ~numel(h.Children)
                    error('No figure exists to add plots too.\nPlease change the ''add2fig'' value to zero.',[]);
                end
                cnt=0;
                for axs = length(h.Children):-1:1
                    axi = h.Children(axs);
                    tmp = whos('axi');                    
                    if strcmp(tmp.class, 'matlab.graphics.axis.Axes')
                        cnt = cnt+1;
                        ax{cnt} = axi;
                    end
                end
                ax2{v} = axes('Position',ax{v}.Position,...
                    'YAxisLocation','right',...
                    'Color','none',...
                    'SortMethod','childorder');
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
        
        if ~strcmp(Plot_,'STI&PESQ')
            arB = area( CI_vec, Res_Bright_area,'LineStyle','none'); hold on;
            %scB = scatter(NoiseLevel,SI_WC_Bright,'bo');
            set(ax2{v},'ColorOrderIndex',ColorInd);
            erB = errorbar(Hrz_Vec,mean(Res_Bright_Matrix),Res_Bright_CI(:,1),Res_Bright_CI(:,2),...
                's' , 'MarkerSize',4); set(erB,'LineWidth',lineWid);
        else
            Res_Bright_CI(:)=0;
            %scB = scatter(NoiseLevel,SI_WC_Bright,'bo');
            set(ax2{v},'ColorOrderIndex',ColorInd);
            erB = scatter(Hrz_Vec,mean(Res_Bright_Matrix),...
                4^2, 's'); set(erB,'LineWidth',lineWid);
        end
            ls='';
            if r_==1, ls='-';
            elseif r_==2, ls='--';
            elseif r_==3, ls=':';
            else ls='-.'; end;
        
            axis([min(Hrz_Vec)-2 20+2 1 4.6*(104/100)]);
        set(ax2{v},'ColorOrderIndex',ColorInd);
        plB = plot(Res_Bright_trend,ls);set(plQ,'LineWidth',lineWid);
        
        %Plot optimal Gdb
        G_=-40:0.1:20;
        pesq_{v}=[Res_Bright_trend(G_)];
        if strcmp(Plot_,'STOI&PESQ'), ls='--';end
        if mod(v,2)==1 %Flat Mask
            SIC = contrast_{v};
            [~,Gdb] = max(SIC);
            op=plot([G_(Gdb) G_(Gdb)],[-200 200],['k' ls]);set(op,'LineWidth',lineWid);
        else % Zone Weighted Mask             
            SIC = contrast_{v}; if strcmp(Plot_,'STOI&PESQ'), SIC = SIC/100;end
            Q = (pesq_{v}-1)/(4.56-1); 
            Lam=1; if strcmp(Plot_,'STI&PESQ'), Lam = 1;end
            [~,Gdb] = max(SIC + Q*Lam);
            op=plot([G_(Gdb) G_(Gdb)],[-200 200],['k' ls]);set(op,'LineWidth',lineWid);
        end
                
        
        if ~add2fig, title(Titles{v}); end
            axis([min(Hrz_Vec)-2 20+2 1 4.6*(104/100)]);
        if ~add2fig, grid on; end;
        if v==length(Version) && strcmp(Plot_,'STOI&PESQ')
            set(gca,'ColorOrderIndex',1);
            stoiB = scatter(0,-10,'o'); %something plotted off the axis
            stoiQ = scatter(0,-10,'r^'); %something plotted off the axis
            set(gca,'ColorOrderIndex',ColorInd);
            pesqB = scatter(0,-10,'s'); %something plotted off the axis
            leg_loc = 'southwest';
            leg = legend([stoiB, stoiQ, pesqB], {'BZ STOI';'QZ STOI';'BZ PESQ'},...
                'Location',leg_loc,...
                'color','w', ...
                'units','pixels',...
                'fontsize',5);
            leg_pos = get(leg,'Position');
            set(leg,'Position', leg_pos + [-37 -62 0 0]);
            
        elseif v==1 && r_==length(rs) && strcmp(Plot_,'STI&PESQ')
                %Measures
                set(ax2{v},'ColorOrderIndex',1);
                stoiB = scatter(0,-10,'o'); %something plotted off the axis
                stoiQ = scatter(0,-10,'r^'); %something plotted off the axis
                set(ax2{v},'ColorOrderIndex',ColorInd);
                pesqB = scatter(0,-10,'s'); %something plotted off the axis
                leg_loc = 'northeast';
                leg = legend([stoiB, stoiQ, pesqB],...
                    {'BZ STI';'QZ STI';'BZ PESQ'},...
                    'Location',leg_loc,...
                    'color','w', ...
                    'units','pixels',...
                    'fontsize',5);
                leg_pos = get(leg,'Position');
                set(leg,'Position', leg_pos + [26 9 0 0]);
        elseif v==3 && r_==length(rs) && strcmp(Plot_,'STI&PESQ')
                %Rooms
                room1 = plot(-100,-10,'k-','LineWidth',lineWid); %something plotted off the axis
                room2 = plot(-100,-10,'k--','LineWidth',lineWid); %something plotted off the axis
                room3 = plot(-100,-10,'k:','LineWidth',lineWid);
                leg_loc = 'northeast';
                leg = legend([room1, room2, room3],...
                    {'Room 1';'Room 2';'Room 3'},...
                    'Location',leg_loc,...
                    'color','w', ...
                    'units','pixels',...
                    'fontsize',5);
                leg_pos = get(leg,'Position');
                set(leg,'Position', leg_pos + [26 -7 0 0]);
        else
            legend off
        end
        
        if ~add2fig, 
            ylim( [1 4.6*(104/100)] );
            ax2{v}.YTick= linspace(1,4.6,6);
            ax2{v}.YTickLabel=num2cell(round(linspace(1,4.6,6),1),1);
        end
           
        if mod(v,2)==0 & r_ == 1
            ylabel({'PESQ (MOS)'});
        else            
            ylabel('');
        end
        
        if ~add2fig, xlabel({'Noise Mask (dB)'}); end; %; '(with reference to Quiet Zone)'});
        
        arB(1).FaceColor = 'none';
        arB(2).FaceColor = plB.Color; 
        drawnow; pause(0.05);  % this is important for transparency!
        arB(2).Face.ColorType = 'truecoloralpha';
        arB(2).Face.ColorData(4) = 0.2*255;              
        
        if add2fig            
            axis([min(Hrz_Vec)-2 20+2 1 4.6*(104/100)]);
            ax2{v}.YTick= linspace(1,4.6,6);
            if r_ == 1 & mod(v,2)==0
                ax2{v}.YTickLabel=num2cell(round(linspace(1,4.6,6),1),1);
            else
                ax2{v}.YTickLabel='';
            end
            ax2{v}.TickLength = [0 0];         
            ax2{v}.XTick=[];
            ax2{v}.XLabel=[];
            ax2{v}.Title=[];
        end
        if mod(v,2)~=0
            set(gca,'YTickLabel',[]);
        end
        
        set(ax2{v},    'FontSize',FontSize-1, ...
            'FontName',Font);

        hold off;
        
    end
    
end
end

if strcmp(Plot_,'STOI&PESQ')
    a= h.Children(2:end);
    leg=h.Children(1);
    h.Children = [a(2:length(a)/2); ...
        a(length(a)/2+2:end); ...
        leg; ...
        a(1); a(length(a)/2+1);];
end

%%
tightfig(h);

set(gcf, 'Units','centimeters', 'Color','w');
fig_pos = get(gcf,'Position');
set(gcf, 'PaperUnits','centimeters', ...
    'Position', [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio], ...
    'PaperSize', [plot_width plot_width/aspect_ratio]);

tightfig(h);

if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
print(['-d' print_fmt], [DocumentPath '\SIC_' Plot_ '_room_' room '_matlab.pdf']);
if strcmp(Plot_,'STI&PESQ'), room='All';end
export_fig([DocumentPath '\SIC_' Plot_ '_room_' room '.pdf'], '-transparent');

%close all;
% Update latex File Name DataBase
Tools.MiKTeX_FNDB_Refresh;

%close all;

