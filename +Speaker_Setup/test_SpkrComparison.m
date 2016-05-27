%clc;
clear;
%close all;
tic;

%%
results_MSEdb = [];
results_contrastdb = [];
results_perc = [];
fprintf('Completed: '); n=0;h=[];
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
%             N_spkrs = 16; 
%             N_spkrs = 24;
            N_spkrs = 32;
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
        [n,h]=Tools.showTimeToCompletion(((spkr-1)*length(freqs)+f)/(length(freqs)*2),n,h);
    end
end
%%
%close all
% [~,ind]=max(results_perc);
% ind1 = logical(mod(ind,2));
% ind2 = [~ind1(2:end) true];
% figure(1);
% plot(freqs,results_perc'); hold on; set(gca,'ColorOrderIndex',1);
% % plot(freqs(ind1),results_perc(1,ind1),'LineWidth',3); set(gca,'ColorOrderIndex',2);
% % plot(freqs(ind2),results_perc(2,ind2),'LineWidth',3); hold off;
% xlim([1e2 1e4]);
% grid on; grid minor
% set(gca,'XScale','log');
% title('Multizone vs Parametic - Percent Difference  (Third Octave)');
% legend({'MSR';'PAL'},'location','southwest');
% ylabel('Percent Level Difference (%)');
% xlabel('Frequency (Hz)');

MSRcol = [0 0.8 0];
PLAcol = [1 0 0];
HYBcol = [0.2 0.2 1];

% figure(str2num(strrep(num2str('Fil'*1),' ','')))
% n = 6; %filter order
% f_c = 1e3;
% G_low   = sqrt(  1./(1+(freqs/f_c).^(2*n))  );G_low = G_low.*G_low; % cascade butterworth
% G_high  = sqrt(  1./(1+1./(freqs/f_c).^(2*n))  );G_high = G_high.*G_high; % cascade butterworth
% G_high2  = sqrt(  1./(1+1./(freqs/(f_c+1e3)).^(2*n))  );G_high2 = G_high2.*G_high2; % cascade butterworth
% G_high3  = sqrt(  1./(1+1./(freqs/(f_c+2e3)).^(2*n))  );G_high3 = G_high3.*G_high3; % cascade butterworth
% G_high4  = sqrt(  1./(1+1./(freqs/(f_c+3e3)).^(2*n))  );G_high4 = G_high4.*G_high4; % cascade butterworth
% plot([f_c f_c],[-1 1]*1e2,'-.k','LineWidth',1,'Color',[0 0 0]); hold on;
% plot([f_c f_c]+1e3,[-1 1]*1e2,'-.k','LineWidth',1,'Color',1-(1-[0 0 0])*0.8);
% plot([f_c f_c]+2e3,[-1 1]*1e2,'-.k','LineWidth',1,'Color',1-(1-[0 0 0])*0.6);
% plot([f_c f_c]+3e3,[-1 1]*1e2,'-.k','LineWidth',1,'Color',1-(1-[0 0 0])*0.4);
% mean1 = mean(mag2db(G_low + G_high));
% mean2 = mean(mag2db(G_low + G_high2));
% mean3 = mean(mag2db(G_low + G_high3));
% mean4 = mean(mag2db(G_low + G_high4));
% plXMEAN4=plot(freqs([1 end]),[mean4 mean4],'--','LineWidth',1,'Color',1-(1-HYBcol)*0.4); 
% plXMEAN3=plot(freqs([1 end]),[mean3 mean3],'--','LineWidth',1,'Color',1-(1-HYBcol)*0.6);
% plXMEAN2=plot(freqs([1 end]),[mean2 mean2],'--','LineWidth',1,'Color',1-(1-HYBcol)*0.8); 
% plXMEAN=plot(freqs([1 end]),[mean1 mean1],'--','LineWidth',1,'Color',HYBcol); 
% plot(freqs,mag2db(G_low + G_high4),'-','LineWidth',1,'Color',1-(1-HYBcol)*0.4);
% plot(freqs,mag2db(G_low + G_high3),'-','LineWidth',1,'Color',1-(1-HYBcol)*0.6); 
% plot(freqs,mag2db(G_low + G_high2),'-','LineWidth',1,'Color',1-(1-HYBcol)*0.8);
% plXSUM=plot(freqs,mag2db(G_low + G_high),'-','LineWidth',1,'Color',HYBcol); 
% plHPF4=plot(freqs,mag2db(G_high4),':','LineWidth',3,'Color',1-(1-PLAcol)*0.4);
% plHPF3=plot(freqs,mag2db(G_high3),':','LineWidth',3,'Color',1-(1-PLAcol)*0.6);
% plHPF2=plot(freqs,mag2db(G_high2),':','LineWidth',3,'Color',1-(1-PLAcol)*0.8);
% plHPF=plot(freqs,mag2db(G_high),':','LineWidth',3,'Color',PLAcol);
% plLPF=plot(freqs,mag2db(G_low),':','LineWidth',3,'Color',MSRcol); 
% hold off; set(gca,'XScale','log');
% xlim([1e2 1e4]); ylim([-70 5]);
% grid on; grid minor;
% set(gca,'fontname','fixedwidth');
% title('12^{th} Order Linkwitz-Riley Crossover Filters','fontname','cambria');
% leg = legend(...
%     [plLPF,plHPF,plHPF2,plHPF3,plHPF4,plXSUM,plXMEAN,plXMEAN2,plXMEAN3,plXMEAN4],...
%     {'LPF@1kHz','HPF@1kHz','HPF@2kHz','HPF@3kHz','HPF@4kHz','LPF+HPF','Mean 1','Mean 2','Mean 3','Mean 4'},'location','southwest','fontname','cambria');
% ylabel('Magnitude (dB)','fontname','cambria');
% xlabel('Frequency (Hz)','fontname','cambria');
% ax=gca; ax.YRuler.MinorTick=ax.YRuler.MinorTick;%force minor ticks to stay as they are


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

[~,ind]=max(results_contrastdb);
ind1 = logical(mod(ind,2));
ind2 = [~ind1(2:end) true];
figure(str2num(strrep(num2str('Con'*1),' ',''))); [ax,h1,h2]=plotyy(1,1,1,1); delete(h2); hold on;
plot([f_cutoff f_cutoff],[-1 1]*1e4,'-.k','LineWidth',1);
% plot([1e2 1e4],[MSR_mean MSR_mean],'--','Color',MSRcol);
% plot([1e2 1e4],[PLA_mean PLA_mean],'--','Color',PLAcol);
% plot([1e2 1e4],[Hybrid_mean Hybrid_mean],'--','Color',HYBcol); 
% plot(freqs,result_filt_low,':','LineWidth',1.5,'Color',MSRcol); 
% plot(freqs,result_filt_high,':','LineWidth',1.5,'Color',PLAcol);
plHYB=plot(freqs,hybrid_contrast,'--','LineWidth',2.5,'Color',[HYBcol 1]);  
plMSR=plot(freqs,results_contrastdb(1,:),'-','LineWidth',2,'Color',[MSRcol 0.5]); 
plPLA=plot(freqs,results_contrastdb(2,:),'-','LineWidth',2,'Color',[PLAcol 0.5]);
hold off;
ax(1).YTickLabelMode = 'manual';
ax(2).YTickLabelMode = 'manual';
ax(1).YTickMode = 'manual';
ax(2).YTickMode = 'manual';
xlim([1e2 1e4]); ylim([0 90]);
ax(1).YColor = [0 0 0];
ax(1).YTick = [0:10:90];
ax(1).YTickLabel = num2str([0:10:90]','%.1f');
grid on; grid minor;
set(ax(1),'XScale','log');
ax(1).FontName = 'fixedwidth';
ax(2).FontName = 'fixedwidth';
title(ax(1),'Third Octave Acoustic Contrast for Reproduction Methods','fontname','cambria');
leg = legend(ax(1),[plMSR,plPLA,plHYB],{'MSR';'PLA';'Hybrid'},'location','northwest','fontname','cambria');
ylabel(ax(1),'Acoustic Contrast (dB)','fontname','cambria');
xlabel(ax(1),'Frequency (Hz)','fontname','cambria');
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
hold off;
ax(1).YRuler.MinorTickValuesMode='manual';%force minor ticks to stay as they are
ax(2).YRuler.MinorTickValuesMode='manual';%force minor ticks to stay as they are

% print('M:\MSR\tex\latex\Parametric\Acoustic_Contrast','-dpdf');

% [~,ind]=max(results_MSEdb);
% ind1 = logical(mod(ind,2));
% ind2 = [~ind1(2:end) true];
% figure(3);
% plot(freqs,results_MSEdb'); hold on; set(gca,'ColorOrderIndex',1);
% % plot(freqs(ind1),results_contrastdb(1,ind1),'LineWidth',3); set(gca,'ColorOrderIndex',2);
% % plot(freqs(ind2),results_contrastdb(2,ind2),'LineWidth',3); hold off;
% xlim([1e2 1e4]);
% grid on; grid minor
% set(gca,'XScale','log');
% title('Multizone vs Parametic - Bright MSE  (Third Octave)');
% legend({'MSR';'PAL'},'location','southwest');
% ylabel('Mean Squared Error (dB)');
% xlabel('Frequency (Hz)');

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script