function SpkrComparison

clc;
clear;
close all;
tic;

%%
%Figure Output Settings
DocumentPath = 'tex\latex\APSIPA2016';
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = 88.9/10 + 6.35/10 + 88.9/10; %IEEE full text width
aspect_ratio = 11/3;
FontSize = 9;
Font = 'Times';
lineWid = 0.5;

%%
axGridSep = [0.5 0.5];
hCon  = figure(str2num(strrep(num2str('Con'*1),' ','')));
hCon_ = tightPlots(1,4,plot_width,[1 1], axGridSep.*[0 1], 1, 1,'centimeters');
hMSE  = figure(str2num(strrep(num2str('MSE'*1),' ','')));
hMSE_ = tightPlots(1,4,plot_width,[1 1], axGridSep.*[0 1], 1, 1,'centimeters');
hALL  = figure(str2num(strrep(num2str('ALL'*1),' ','')));
hALL_ = tightPlots(2,4,plot_width,[1 1], axGridSep.*[1 1], 1, 1,'centimeters');

%%
results_MSEdb = [];
results_contrastdb = [];
results_perc = [];
leg=matlab.graphics.GraphicsPlaceholder;
fprintf('\tCompleted: '); n=0;

maxL = ceil(pi*(2*(2*pi*8000/343*0.9)+1)/2/pi)+1;
N_spkrs_ = [16, 24, 32, maxL];
for nspkr = 1:length(N_spkrs_)
    N_spkrs = N_spkrs_(nspkr);
    for spkr = 1:2
        %freqs = (10.0^3) * ((2.0) .^ ([-8:9]./3)); %Third octave
        freqs = logspace(log10(3000),log10(8000),100);
        for f = 1:length(freqs)
            
            %%
            room_width = 3.200; % metres
            rad = ( room_width - 2*0.115 - 2*0.085 ) / 2; % metres
            
            % speech_layout = {'brightzone_pos_angle',        180, ...
            %                  'quietzone_pos_angle',         0, ...
            %                  'brightzone_source_angle',     14.5};
            % masker_layout = {'brightzone_pos_angle',        0, ...
            %                  'quietzone_pos_angle',         180, ...
            %                  'brightzone_source_angle',     0};
            
            speech_layout = {'brightzone_pos_angle',        90, ...
                'quietzone_pos_angle',         -90, ...
                'brightzone_source_angle',     0};
            masker_layout = { ...
                'brightzone_pos_angle',        -90, ...
                'quietzone_pos_angle',         90, ...
                'brightzone_source_angle',     180, ... -atand(0.6/1.3), ...
                'brightzone_source_type',      'ps'};
            
            array_type = 'circle';
            spkr_radius = 1.3;
            
            if strcmp(array_type, 'circle')
                [x,y] = pol2cart(-90/180*pi, 0.6);
                x_ = sqrt(spkr_radius^2-y^2);
                % %             For a different source angle
                %             x_ = sqrt(1.3^2+0.6^2);
                spkr_spacing = []; %Auto-calculate spacing
            elseif strcmp(array_type, 'line')
                x_=spkr_radius;
                % %             For a different source angle
                %             x_ = sqrt(1.3^2+0.6^2);
                %             th = pi-atan(0.6/1.3);
                spkr_spacing = 0.001; %1mm spacing between adjacent loudspeakers
            end
            Para_Spkr = Parametric_Synthesis.parametric_soundfield;
            if spkr ==1
                masker = { ...
                    'angleto_firstloudspeaker',      90, ...
                    'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
                    'numberof_loudspeakers',         N_spkrs, ...
                    'loudspeaker_model',             'Genelec 8010A', ...
                    'loudspeaker_radius',            spkr_radius, ...
                    'loudspeaker_spacing',           spkr_spacing, ...
                    'speaker_array_type',            array_type, ...
                    'brightzone_source_dist',        x_, ...
                    'angleof_loudspeakerarrcentre', 180, ...
                    'quiet_weight',                 1e2};
            elseif spkr == 2
                %parametric
                masker = {  ...
                    'angleto_firstloudspeaker',     atan2d(-0.6,-x_), ...
                    'angleof_loudspeakerarrcentre', 180, ... +atand(0.6/1.3), ...
                    'loudspeaker_radius',           x_, ... spkr_radius, ...
                    'numberof_loudspeakers',        1, ...
                    'loudspeaker_model',            'Parametric', ...
                    'loudspeaker_spacing',          0.01, ...
                    'speaker_array_type',           'line', ...
                    'brightzone_source_dist',        x_};
                Para_Spkr = Para_Spkr.set_f1( 40000 );
            end
            setup = Speaker_Setup.createSetup({...
                'frequency',                    freqs(f), ...
                masker_layout{:}, ...
                masker{:}, ...
                'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
                'reproduction_radius',          1.0, ...
                'bright_weight',                1.0, ...
                'unattended_weight',            0.05, ...
                'brightzone_radius',            0.3, ...
                'brightzone_pos_distance',      0.6, ...
                'quietzone_radius',             0.3, ...
                'quietzone_pos_distance',       0.6, ...
                'maximum_frequency',            8000, ...
                'loudspeaker_object',           Para_Spkr });
            
            %%
            setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
            
            setup = setup.calc_Loudspeaker_Weights();
            setup = setup.reproduceSoundfield('SAMPLES_ONLY');
            
            if spkr == 1
                f_cutoff = Broadband_Tools.getAliasingFrequency(setup, 'old')/2/pi*343;
            end
            
            results_perc(spkr,f)       = (setup.Bright_Sample - setup.Quiet_Sample)*100;
            results_contrastdb(spkr,f) = mag2db(setup.Acoustic_Contrast);
            results_MSEdb(spkr,f) = mag2db(setup.MSE_Bright);
            %%
            [n]=Tools.showTimeToCompletion(((nspkr-1)*2*length(freqs) + (spkr-1)*length(freqs)+f)/(length(freqs)*2*length(N_spkrs_)),n);
        end
    end
    LRorder = 12;
    leg_ = plotComparisons( N_spkrs, nspkr, LRorder, freqs, f_cutoff, results_contrastdb, results_MSEdb, hCon_, hMSE_, hALL_ );
    if nspkr ~= length(N_spkrs_), fprintf(repmat('\n',1,n)); end
    if ~isempty(leg_) && ~isa(leg_,'matlab.graphics.GraphicsPlaceholder')
        if length(leg)==1 && isa(leg,'matlab.graphics.GraphicsPlaceholder')
            leg=leg_;
        else
            leg(end+1:end+length(leg_))=leg_;
        end
    end
end



%% Separated Figures
saveComparisons(hCon, leg(1), ...
    'Acoustic Contrast for Reproduction Methods', ...
    DocumentPath, plot_width, aspect_ratio, [0.02 +0.05], -0.08, 1)

saveComparisons(hMSE, leg(3), ...
    'Mean Squared Error for Reproduction Methods', ...
    DocumentPath, plot_width, aspect_ratio, [0.04 +0.05], -0.08, 1)

%% Combined Figures
for s = 1:length(leg(2).String)
leg(2).String{s} = strrep(leg(2).String{s},'\zeta_','');
end
axs = findobj(hALL.Children,'type','axes'); axPos=[];keepflag=true;
for a = 1:length(axs)
    if ~strcmpi(axs(a).XLabel.String,' ') && ~isempty(axs(a).XLabel.String)
        tmpUnit = axs(a).Units;
        axs(a).Units = 'Centimeters';
        axPos(end+1,:) = axs(a).Position(1);
        axs(a).Units = tmpUnit;
        if ~keepflag
            axs(a).XLabel.String = ' ';
        else
            axs_no = a; 
        end
        keepflag = false;
    end
end
tmpUnit1 = axs(axs_no).XLabel.Units; tmpUnit2 = axs(axs_no).Units;
axs(axs_no).XLabel.Units = 'Centimeters'; axs(axs_no).Units = 'Centimeters';
axs(axs_no).XLabel.Position(1) = axs(axs_no).XLabel.Position(1) + mean(axPos,1) ...
    -axs(axs_no).Position(1);
axs(axs_no).XLabel.Units = tmpUnit1; axs(axs_no).Units = tmpUnit2;
leg(4).Visible = 'off';

saveComparisons(hALL, leg([2 4]), ...
    'Acoustic Contrast and Mean Squared Error for Reproduction Methods', ...
    DocumentPath, plot_width, aspect_ratio, [0.02 -0.02; 0.04 +0.01], -0.06, 2)


%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

end

function  saveComparisons(hFig, hLeg, figTitle, DocumentPath, plot_width, aspect_ratio, legshift, titleshift, rows)
figure(hFig)

tightfig(hFig);
hFig.Units = 'centimeters';
hFig.Color = 'w';
fig_pos = hFig.Position;
hFig.Position = [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio*rows];pause(0.05);drawnow;
hFig.PaperUnits = 'centimeters';
hFig.PaperSize = [plot_width plot_width/aspect_ratio*rows];
tightfig(hFig);
hFig.Position = hFig.Position + [0 0 0 1];pause(0.05);drawnow;
for l = 1:length(hLeg)
hLeg(l).Position = hLeg(l).Position + [legshift(l,1) legshift(l,2) 0 0];
end

ht=suptitle(figTitle);
ht.FontName = 'times';
ht.FontSize = 10;
ht.Position = ht.Position + [0 titleshift 0];


if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
export_fig([DocumentPath '\' figTitle '.pdf']);
end


function [leg] = plotComparisons( N_spkrs, nspkr, LRorder, freqs, f_cutoff, results_contrastdb, results_MSEdb, hCon_, hMSE_, hALL_ )

%%
MSRcol = [0 0.8 0];
PLAcol = [1 0 0];
HYBcol = [0.2 0.2 1];
% MSRcol = [77 175 74]/255;
% PLAcol = [228 26 28]/255;
% HYBcol = [55 126 184]/255;

G_low   = sqrt(  1./(1+(freqs/f_cutoff).^LRorder)  );
G_low = G_low.*G_low; % cascade butterworth
G_high  = sqrt(  1./(1+1./(freqs/f_cutoff).^LRorder)  );
G_high = G_high.*G_high; % cascade butterworth

% Acoustic Contrast
result_filt_low  = mag2db(db2mag(results_contrastdb(1,:)/2).*G_low)*2;
result_filt_high = mag2db(db2mag(results_contrastdb(2,:)/2).*G_high)*2;
hybrid_contrast = mag2db(db2mag(result_filt_low/2)+db2mag(result_filt_high/2))*2;

% Mean Squared Error
result_filt_low  = mag2db(db2mag(results_MSEdb(1,:)).*G_low);
result_filt_high = mag2db(db2mag(results_MSEdb(2,:)).*G_high);
hybrid_MSE = mag2db(db2mag(result_filt_low)+db2mag(result_filt_high));

MSR_AC_mean = mean(results_contrastdb(1,:));
PLA_AC_mean = mean(results_contrastdb(2,:));
Hybrid_AC_mean = mean(hybrid_contrast);
MSR_MSE_mean = mean(results_MSEdb(1,:));
PLA_MSE_mean = mean(results_MSEdb(2,:));
Hybrid_MSE_mean = mean(hybrid_MSE);

fprintf('\n');
fprintf('\tNumber of Speakers:   %d\n',     N_spkrs);
fprintf('\tMSR Mean Contrast:    %.1fdB\n', MSR_AC_mean);
fprintf('\tPL Mean Contrast:     %.1fdB\n', PLA_AC_mean);
fprintf('\tHybrid Mean Contrast: %.1fdB\n', Hybrid_AC_mean);
fprintf('\tMSR Mean MSE:         %.1fdB\n', MSR_MSE_mean);
fprintf('\tPL Mean MSE:          %.1fdB\n', PLA_MSE_mean);
fprintf('\tHybrid Mean MSE:      %.1fdB\n', Hybrid_MSE_mean);

leg = matlab.graphics.GraphicsPlaceholder;

leg_ = createPlot( N_spkrs, freqs, f_cutoff, ...
    results_contrastdb, hybrid_contrast, ...
    [MSRcol;PLAcol;HYBcol], hCon_, nspkr, 'Acoustic Contrast', '\zeta', 1, 'A', true);
if ~isempty(leg_),
    if length(leg)==1 && isa(leg,'matlab.graphics.GraphicsPlaceholder')
        leg=leg_;
    else
        leg(end+1)=leg_;
    end
end

leg_ = createPlot( N_spkrs, freqs, f_cutoff, ...
    results_MSEdb, hybrid_MSE, ...
    [MSRcol;PLAcol;HYBcol], hMSE_, nspkr, 'Mean Squared Error', '\epsilon', 4, 'A', true);
if ~isempty(leg_),
    if length(leg)==1 && isa(leg,'matlab.graphics.GraphicsPlaceholder')
        leg=leg_;
    else
        leg(end+1)=leg_;
    end
end

leg_ = createPlot( N_spkrs, freqs, f_cutoff, ...
    results_contrastdb, hybrid_contrast, ...
    [MSRcol;PLAcol;HYBcol], hALL_(1:4), nspkr, 'Acoustic Contrast', '\zeta', 1, 'A', false);
if ~isempty(leg_),
    if length(leg)==1 && isa(leg,'matlab.graphics.GraphicsPlaceholder')
        leg=leg_;
    else
        leg(end+1)=leg_;
    end
end

leg_ = createPlot( N_spkrs, freqs, f_cutoff, ...
    results_MSEdb, hybrid_MSE, ...
    [MSRcol;PLAcol;HYBcol], hALL_(5:end), nspkr, 'Mean Squared Error', '\epsilon', 4, 'A'+4, true);
if ~isempty(leg_),
    if length(leg)==1 && isa(leg,'matlab.graphics.GraphicsPlaceholder')
        leg=leg_;
    else
        leg(end+1)=leg_;
    end
end

end

function leg = createPlot( N_spkrs, freqs, f_cutoff, results, hybrid_result, colors, axes_handle, axes_no, yaxeslbl, measSym, leg_axis, LetterStart, showXlbl)
MSRcol = colors(1,:);
PLAcol = colors(2,:);
HYBcol = colors(3,:);
[~,ind]=max(results);
ind1 = logical(mod(ind,2));
ind2 = [~ind1(2:end) true];
%     subplot(1,4,nspkr);
axes(axes_handle(axes_no));
[ax,h1,h2]=plotyy(1,1,1,1); delete(h2); hold on;
plCUT = plot([f_cutoff f_cutoff],[-1 1]*1e4,'-.k','LineWidth',1);
plHYB = plot(freqs,hybrid_result,'-','LineWidth',1.5,'Color',[HYBcol 1]);
plMSR = plot(freqs,results(1,:),'--','LineWidth',1.0,'Color',[MSRcol 1]);
plPLA = plot(freqs,results(2,:),'--','LineWidth',1.0,'Color',[PLAcol 1]);
hold off;
ax(1).TickDirMode = 'manual';
ax(1).TickDir = 'both';
ax(1).YTickLabelMode = 'manual';
ax(2).YTickLabelMode = 'manual';
ax(1).YTickMode = 'manual';
ax(2).YTickMode = 'manual';
xlim([1e2 20^3]);
ax(1).XTick = [10^2 10^3 20^3];
ax(1).XTickLabelMode = 'manual';
if showXlbl
    ax(1).XTickLabel = {'0.1' '1' '8'};
else
    ax(1).XTickLabel = {};
end
ax(1).TickLength = ax(1).TickLength .* [1.5 1];
ax(1).YColor = [0 0 0];
switch yaxeslbl
    case 'Acoustic Contrast'        
        limit_vec = [0 150];
        ylim(limit_vec);
        ax(1).YTick = limit_vec(1):20:limit_vec(end);
        if axes_no == 1
            ylbls = mat2cell(num2str(ax(1).YTick','%.0f'), ones(length(ax(1).YTick),1), 3);
            for yl = 1:length(ylbls)
                ylbls{yl} =  strrep(ylbls{yl},' ','');
            end
            ax(1).YTickLabel = ylbls;
        else
            ax(1).YTickLabel = '';
        end
        text(120,limit_vec(end)-0.1*diff(limit_vec),['(' char(64+axes_no) ')'],'fontname','times','fontsize',10)
        text(120,limit_vec(end)-0.25*diff(limit_vec),['$L = ' num2str(N_spkrs) '$'],'interpreter','latex','fontname','times','fontsize',10)
        
    case 'Mean Squared Error'
        limit_vec = [-50 0];
        y_inc = 10;
        ylim(limit_vec);
        ax(1).YTick = limit_vec(1):y_inc:limit_vec(end);
        if axes_no == 1
            ylbls = mat2cell(num2str(ax(1).YTick','%.0f'), ones(length(ax(1).YTick),1), 3);
            for yl = 1:length(ylbls)
                ylbls{yl} =  strrep(ylbls{yl},' ','');
            end
            ax(1).YTickLabel = ylbls;
        else
            ax(1).YTickLabel = '';
        end
        text(120,limit_vec(end)-0.1*diff(limit_vec),['(' char(LetterStart-1+axes_no) ')'],'fontname','times','fontsize',10)
        text(120,limit_vec(end)-0.25*diff(limit_vec),['$L = ' num2str(N_spkrs) '$'],'interpreter','latex','fontname','times','fontsize',10)
end
grid on; grid minor;
set(ax(1),'XScale','log');
ax(1).FontName = 'fixedwidth';
ax(1).FontSize = 10;
ax(2).FontName = 'fixedwidth';
ax(2).FontSize = 10;
if axes_no == 1
    ylabel(ax(1),[yaxeslbl ' (dB)'],'fontname','times','fontsize',10);
end
if showXlbl
    xlabel(ax(1),'Frequency (kHz)','fontname','times','fontsize',10);
end
leg=[];
if axes_no == leg_axis
    [leg, legIcons] = legend(ax(1),[plMSR,plPLA,plHYB,plCUT],...
        {['$' measSym '_{\mathrm{MSR}}$'];['$' measSym '_{\mathrm{PL}}$'];['$' measSym '_{\mathcal{H}}$'];'$k_{\mathrm{u}}$'}, ...
        'box','off','color','none', ...
        'location','northeast','fontname','times','fontsize',9,'interpreter','latex');
    leglines = findobj(legIcons,'Type','line');
    for ll = 1:2:length(leglines)
        leglines(ll).XData(1)  = leglines(ll).XData(1) + diff(leglines(ll).XData) .* 5/12;
    end
end
hold on;
ax(2).YAxisLocation = 'left';
ax(2).Color = 'none';
ax(2).YColor = ax(1).YColor;
ax(2).Box = 'off';
ax(2).TickLength = [0 0];
ax(2).XTick = [];
% ax(2).YTick = sort(round([MSR_mean+1 PLA_mean-1 Hybrid_mean],0));
ax(2).YLim = [0 90];
% ax(2).YTickLabel = num2str(sort([MSR_mean PLA_mean Hybrid_mean]'),'%.1f');
ax(2).YTick = [];
ax(2).YTickLabel = '';
hold off;
%     ax(1).YRuler.MinorTickValuesMode='manual';%force minor ticks to stay as they are
%     ax(2).YRuler.MinorTickValuesMode='manual';%force minor ticks to stay as they are
end