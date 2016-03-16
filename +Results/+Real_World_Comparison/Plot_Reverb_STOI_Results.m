clc;
clear;
close all;

%%
ColorInd = 4; %this is for whatever the PESQ colour will be plotted as

%% Info
results_type = 'STOI';

Drive = 'Z:\'; % Database drive (storage drive)

%LUT_resolution =  '512f_256w'; %Look-Up Table resolution
LUT_resolution =  '512f_32w'; %Look-Up Table resolution
signal_info.c = 343; % Speed of sound in metres/sec
signal_info.Fs = 16000; % Sampling frequency
signal_info.Nfft = 1024;% Number of fft components
signal_info.overlap = 0.5;
signal_info.f_low  = 150;  % Hz
signal_info.f_high = 8000; % Hz
signal_info.weight = [];
signal_info.L_noise_mask = []; % dB
signal_info.input_filename = [];

array_type = 'circle';
spkr_radius = 1.3;

parametric_speaker = Parametric_Synthesis.parametric_soundfield;
parametric_speaker.P1 = db2mag( 100 ); % 100dB amplitude parametric array loudspeaker
parametric_speaker.P2 = db2mag( 100 ); % 100dB secondary amplitude
if strcmp(array_type, 'circle')
    [x,y] = pol2cart(-90/180*pi, 0.6);
    x_ = sqrt(spkr_radius^2-y^2);
    th_c = atan2(y,-x_);
    th = th_c;
    spkr_spacing = []; %Auto-calculate spacing
elseif strcmp(array_type, 'line')
    x_=spkr_radius;
    th_c = 180;
    th = atan2(-0.6,-spkr_radius);
    spkr_spacing = 0.001; %1mm spacing between adjacent loudspeakers
end
N_spkrs = 24;
loudspeaker_layout = { ...
    'angleto_firstloudspeaker',      90, ...
    'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
    'numberof_loudspeakers',         N_spkrs, ...
    'loudspeaker_model',             'Genelec 8010A', ...
    'loudspeaker_radius',            spkr_radius, ...
    'loudspeaker_spacing',           spkr_spacing, ...
    'speaker_array_type',            array_type };
other = { ...
    'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
    'reproduction_radius',          1.0, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 0, ...
    'unattended_weight',            0.05, ...
    'brightzone_radius',            0.3, ...
    'brightzone_pos_distance',      0.6, ...
    'quietzone_radius',             0.3, ...
    'quietzone_pos_distance',       0.6, ...
    'maximum_frequency',            8000, ...
    'angleof_loudspeakerarrcentre', 180, ...
    'loudspeaker_object',           parametric_speaker };
speech_layout = { ...
    'brightzone_pos_angle',        90, ...
    'quietzone_pos_angle',         -90, ...
    'brightzone_source_angle',     0, ...
    'brightzone_source_type',      'pw'};
speech_loudspeaker_layout = { ...
    'angleto_firstloudspeaker',      90, ...
    'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
    'numberof_loudspeakers',         N_spkrs, ...
    'loudspeaker_model',             'Genelec 8010A', ...
    'loudspeaker_radius',            spkr_radius, ...
    'loudspeaker_spacing',           spkr_spacing, ...
    'speaker_array_type',            array_type, ...
    'brightzone_source_dist',        x_, ...
    other{:}};

Setup = Speaker_Setup.createSetup({ speech_layout{:}, speech_loudspeaker_layout{:}});

Room_Setup = Room_Acoustics.Room;
Room_Setup.NoReceivers = 32;





% LUT_resolution = '512f_256w';
% loudspeakers   = 295;  % Number of loudspeakers
% speaker_arc    = 360;  % Degrees
% speaker_radius = 1.5; % Metres
% Num_Receivers = 32;

% pw_angle = {0; ...
%     0; ...
%     15; ...
%     15; ...
%     90; ...
%     90;};

mask_type_ = {'Zone Weight Masker AliasCtrl'; ...
    'Zone Weight Masker AliasCtrl';};

for i=1:length(mask_type_)
    Version{1,i} = strrep(strrep(mask_type_{i},'Weighted','Weight'),' ','');
end

h=figure(1);
subp = [];
rs = [ 1 ];

for r_ = 1:length(rs)
    r = rs(r_);
    if r == 1
        % % ROOM 1
        % % Anechoic6
        Room_Setup.setRoomSize( [10 10 10] ); %Anechoic
        Room_Setup.Wall_Absorb_Coeff = 1.0;
    elseif r == 2
        % % ROOM 2
        % % Small Open Plan Office
        Room_Setup.setRoomSize( [4 9 3] );   %Small Open Plan Office
        Room_Setup.Wall_Absorb_Coeff = 0.3;
    elseif r == 3
        % % ROOM 3
        % % Medium Open Plan Office
        Room_Setup.setRoomSize( [8 10 3] );   %Medium Open Plan Office
        Room_Setup.Wall_Absorb_Coeff = 0.3;
    elseif r == 4
        % % ROOM 4
        % % Cafe / Restaurant
        Room_Setup.setRoomSize( [9 14 3] );   %Cafe/Restaurant
        Room_Setup.Wall_Absorb_Coeff = 0.3;
    end
    
    Room_Setup.setReproductionCentre( Room_Setup.Room_Size / 2 ); %[5 5 5];
    if Room_Setup.Wall_Absorb_Coeff >= 1.0
        room_type = 'Anechoic';
    else
        room_type = 'Reverberant';
    end
    
    %     Version = {['10000weight__with' mask_type{1}]; ...
    %         ['10000weight__with' mask_type{2}]};
    
    Titles  = { ...
        ['' mask_type_{1} ]; ... ' and \theta=' num2str(pw_angle{1}) '°']; ...
        ['' mask_type_{2} ]};% ' and \theta=' num2str(pw_angle{6}) '°']};
    
    
    %Figure Output Settings
    DocumentPath = 'tex\latex\Parametric';
    print_fmt = 'pdf'; %figure image file format
    print_res = 600; %dpi
    plot_width = (88.9)/10;% + 6.35 + 88.9)/10; %IEEE full text width
    aspect_ratio = 3/3;
    FontSize = 9;
    Font = 'Times';
    lineWid = 0.5;
    
    
    %     %%
    %     room = strrep(sprintf(strrep(repmat('%g',1,length(Room_Size)),'g%','g %'),Room_Size),' ','x');
    %     room_cent = strrep(sprintf(strrep(repmat('%g',1,length(Reproduction_Centre)),'g%','g %'),Reproduction_Centre),' ','x');
    
    %%
    
    %h = figure(1);
    for v = 1:length(Version)
        %         %% Create paths
        %         ResultsPath = ['+Results\+Reverb__' num2str(Num_Receivers) 'Rec_' room 'Dim_' room_cent 'Ctr_' num2str(Wall_Absorption_Coeff) 'Ab\'];
        %         Output_file_path_ext = ['+' num2str(speaker_radius*2) 'm_SpkrDia\+' num2str(loudspeakers) 'Spkrs_' num2str(speaker_arc) 'DegArc_LUT_' LUT_resolution '\'];
        %         Results_filepath = [ResultsPath Output_file_path_ext 'STI_Results_' num2str(pw_angle{v}) 'pwAngle_' Version{v} '.csv'];
        
        signal_info.method = Version{v};
        if v==1
            signal_info.recording_type = 'realworld';
            Results_filepath = [ ...
                Results.getResultsPath( Setup, LUT_resolution, Room_Setup, signal_info, Drive ), ...
                results_type '_Results.csv'];
        elseif v==2
            signal_info.recording_type = 'realworld';
            Results_filepath = [ ...
                Results.getResultsPath( Setup, LUT_resolution, Room_Setup, signal_info, Drive ), ...
                results_type '_Results.csv'];
        end
        % Read results
        [NoiseLevel,SI_WC_Bright,SI_WC_Quiet,ConfInt_Bright_Low,ConfInt_Bright_Up,ConfInt_Quiet_Low,ConfInt_Quiet_Up] = Results.import_SpeechIntelligibility_Reverb(Results_filepath);
        
        
        % Generate Plottable data matrices and vectors
        [Hrz_Vec, Res_Bright_Matrix, Res_Bright_trend, Res_Bright_area, Res_Bright_CI, CI_vec] = Results.generatePlotData( NoiseLevel, SI_WC_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up);
        [~, Res_Quiet_Matrix, Res_Quiet_trend, Res_Quiet_area, Res_Quiet_CI, ~] = Results.generatePlotData( NoiseLevel, SI_WC_Quiet, ConfInt_Quiet_Low, ConfInt_Quiet_Up);
        
        
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
            %            text(txt_pos(1),txt_pos(2),['\theta = ' num2str(pw_angle{v}) '°'],'fontname',Font,'fontsize',FontSize)
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
    
end
set(gcf, 'Units','centimeters', 'Color','w');
fig_pos = get(gcf,'Position');
set(gcf, 'PaperUnits','centimeters', ...
    'Position', [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio], ...
    'PaperSize', [plot_width plot_width/aspect_ratio]);
%%
% tightfig(h);
% set(gcf, 'Units','centimeters', 'Color','w');
% fig_pos = get(gcf,'Position');
% set(gcf, 'PaperUnits','centimeters', ...
%     'Position', [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio], ...
%     'PaperSize', [plot_width plot_width/aspect_ratio]);
%
% tightfig(h);
%
% if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
% %print(['-d' print_fmt], [DocumentPath '\SIC_Enhancement']);
% export_fig([DocumentPath '\SIC_STI_AllRooms.pdf'], '-transparent');
%
% %close all;
% % Update latex File Name DataBase
% Tools.MiKTeX_FNDB_Refresh;

%close all;
