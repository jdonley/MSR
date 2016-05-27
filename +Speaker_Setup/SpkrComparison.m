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
h  = figure(str2num(strrep(num2str('Con'*1),' ','')));
ha = tightPlots(1,4,plot_width,[1 1], [0 0.5], 1, 1,'centimeters');

%%
results_MSEdb = [];
results_contrastdb = [];
results_perc = [];
fprintf('\tCompleted: '); n=0;

maxL = ceil(pi*(2*(2*pi*8000/343*0.9)+1)/2/pi)+1;
N_spkrs_ = [16, 24, 32, maxL];
for nspkr = 1:length(N_spkrs_)
    N_spkrs = N_spkrs_(nspkr);
    for spkr = 1:2
        %freqs = (10.0^3) * ((2.0) .^ ([-8:9]./3)); %Third octave
        freqs = logspace(log10(100),log10(8000),100);
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
    %%
    MSRcol = [0 0.8 0];
    PLAcol = [1 0 0];
    HYBcol = [0.2 0.2 1];
    
    G_low   = sqrt(  1./(1+(freqs/f_cutoff).^(2*n))  );
    G_low = G_low.*G_low; % cascade butterworth
    G_high  = sqrt(  1./(1+1./(freqs/f_cutoff).^(2*n))  );
    G_high = G_high.*G_high; % cascade butterworth
    
    result_filt_low  = mag2db(db2mag(results_contrastdb(1,:)/2).*G_low)*2;
    result_filt_high = mag2db(db2mag(results_contrastdb(2,:)/2).*G_high)*2;
    hybrid_contrast = mag2db(db2mag(result_filt_low/2)+db2mag(result_filt_high/2))*2;
    MSR_mean = mean(results_contrastdb(1,:));
    PLA_mean = mean(results_contrastdb(2,:));
    Hybrid_mean = mean(hybrid_contrast);
    
    fprintf('\n');
    fprintf('\tNumber of Speakers:   %d\n',     N_spkrs);
    fprintf('\tMSR Mean Contrast:    %.1fdB\n', MSR_mean);
    fprintf('\tPL Mean Contrast:     %.1fdB\n', PLA_mean);
    fprintf('\tHybrid Mean Contrast: %.1fdB\n', Hybrid_mean);
    fprintf('\tImprovement over MSR: %.1fdB\n', Hybrid_mean-MSR_mean);
    fprintf('\tImprovement over PL:  %.1fdB\n', Hybrid_mean-PLA_mean);
    fprintf(['\n' repmat(' ',1,n)]);
    
    [~,ind]=max(results_contrastdb);
    ind1 = logical(mod(ind,2));
    ind2 = [~ind1(2:end) true];
    %     subplot(1,4,nspkr);
    axes(ha(nspkr));
    [ax,h1,h2]=plotyy(1,1,1,1); delete(h2); hold on;
    plCUT = plot([f_cutoff f_cutoff],[-1 1]*1e4,'-.k','LineWidth',1);
    plMSR = plot(freqs,results_contrastdb(1,:),'-','LineWidth',1.0,'Color',[MSRcol 1]);
    plPLA = plot(freqs,results_contrastdb(2,:),'-','LineWidth',1.0,'Color',[PLAcol 1]);
    plHYB = plot(freqs,hybrid_contrast,'--','LineWidth',1.5,'Color',[HYBcol 1]);
    hold off;
    ax(1).TickDirMode = 'manual';
    ax(1).TickDir = 'both';
    ax(1).YTickLabelMode = 'manual';
    ax(2).YTickLabelMode = 'manual';
    ax(1).YTickMode = 'manual';
    ax(2).YTickMode = 'manual';
    xlim([1e2 1e4]);
    ax(1).XTick = [10^2 10^3 20^3];
    ax(1).XTickLabelMode = 'manual';
    ax(1).XTickLabel = {'0.1' '1' '8'};
    ax(1).TickLength = ax(1).TickLength .* [1.5 1];
    ylim([0 150]);
    ax(1).YColor = [0 0 0];
    ax(1).YTick = [0:20:140];
    if nspkr == 1
        ylbls = mat2cell(num2str([0:20:140]','%.0f'), ones(length([0:20:140]),1), 3);
        for yl = 1:length(ylbls)
            ylbls{yl} =  strrep(ylbls{yl},' ','');
        end
        ax(1).YTickLabel = ylbls;
        %     elseif nspkr == 4
        %         ylim([0 120]);
        %         ax(1).YTick = [0:20:120];
        %         ax(1).YTickLabel = num2str([0:20:120]','%.0f');
    else
        ax(1).YTickLabel = '';
    end
    grid on; grid minor;
    set(ax(1),'XScale','log');
    ax(1).FontName = 'fixedwidth';
    ax(1).FontSize = 10;
    ax(2).FontName = 'fixedwidth';
    ax(2).FontSize = 10;
    if nspkr == 1
        ylabel(ax(1),'Acoustic Contrast (dB)','fontname','times','fontsize',10);
    end
    xlabel(ax(1),'Frequency (kHz)','fontname','times','fontsize',10);
    if nspkr == 1
        [leg, legIcons] = legend(ax(1),[plMSR,plPLA,plHYB,plCUT],...
            {'$\zeta_{\mathrm{MSR}}$';'$\zeta_{\mathrm{PL}}$';'$\zeta_{\mathcal{H}}$';'$k_{\mathrm{u}}$'}, ...
            'box','off','color','none', ...
            'location','northeast','fontname','times','fontsize',9,'interpreter','latex');
        leglines = findobj(legIcons,'Type','line');
        for ll = 1:2:length(leglines)
            leglines(ll).XData(1)  = leglines(ll).XData(1) + diff(leglines(ll).XData) .* 1/3;
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
    text(120,140,['(' char(64+nspkr) ')'],'fontname','times','fontsize',10)
    text(120,110,['$L = ' num2str(N_spkrs) '$'],'interpreter','latex','fontname','times','fontsize',10)
    hold off;
    %     ax(1).YRuler.MinorTickValuesMode='manual';%force minor ticks to stay as they are
    %     ax(2).YRuler.MinorTickValuesMode='manual';%force minor ticks to stay as they are
    
end



%%
tightfig(h);
h.Units = 'centimeters';
h.Color = 'w';
fig_pos = h.Position;
h.Position = [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio];
h.PaperUnits = 'centimeters';
h.PaperSize = [plot_width plot_width/aspect_ratio];
tightfig(h);
h.Position = h.Position + [0 0 0 1];
% leg.Position = leg.Position + [0.01 -0.06 0 0];
leg.Position = leg.Position + [0.02 +0.05 0 0];

ht=suptitle('Acoustic Contrast for Reproduction Methods');
ht.FontName = 'times';
ht.FontSize = 10;
ht.Position = ht.Position + [0 -0.08 0];

if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
% print(['-d' print_fmt], [DocumentPath '\MSR_PL_HYB_AcousticContrast_matlab.pdf']);
export_fig([DocumentPath '\MSR_PL_HYB_AcousticContrast.pdf']);

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script