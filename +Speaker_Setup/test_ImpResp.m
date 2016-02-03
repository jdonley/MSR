%clc;
clear;
%close all;
tic;

%%
results_contrastdb = [];
results_perc = [];
for spkr = 2
    Nfft = 1024;
    Fs = 16000;
    len = Nfft;
    freqs = linspace(0, Fs/2, len/2 + 1);
    freqs = freqs(1:end-1);
    %freqs = logspace(log10(150),log10(8000),100);
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
        masker_layout = {'brightzone_pos_angle',        -90, ...
            'quietzone_pos_angle',         90, ...
            'brightzone_source_angle',     0};
        
        array_type = 'circle';
        if spkr ==1
            masker = {'angleto_firstloudspeaker',      90, ...
                'numberof_loudspeakers',         24, ...
                'loudspeaker_model',             'Genelec 8010A', ... 
                  'loudspeaker_radius',            1.5, ...                 
                  'loudspeaker_spacing',           0.01, ...
                  'speaker_array_type',            array_type};
        elseif spkr == 2
            %parametric
            if strcmp(array_type, 'circle')
                [x,y] = pol2cart(-90/180*pi, 0.6);
                x_ = sqrt(1.5^2-y^2);
                th_c = atan2(y,-x_);
                th = th_c;
            elseif strcmp(array_type, 'line')
                x_=1.5;
                th_c = 180;
                th = atan2(-0.6,-1.5);
            end
            masker = {'angleto_firstloudspeaker',     th/pi*180, ...
                     'loudspeaker_radius',           x_, ...
                     'numberof_loudspeakers',        1, ...
                     'loudspeaker_model',            'Parametric', ...                  
                     'loudspeaker_spacing',           0.01, ...
                     'speaker_array_type',            'line'};
        end
        setup = Speaker_Setup.createSetup({...
            'frequency',                    freqs(f), ...
            masker_layout{:}, ...
            masker{:}, ...
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
            'angleof_loudspeakerarc',       180 , ...
            'loudspeaker_object',           Parametric_Synthesis.parametric_soundfield});
        
        %%
        setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
        
        setup = setup.calc_Loudspeaker_Weights();
        setup = setup.reproduceSoundfield('SAMPLES_ONLY');                
        
        results_perc(spkr,f)       = (setup.Bright_Sample - setup.Quiet_Sample)*100;
        results_contrastdb(spkr,f) = mag2db(setup.Acoustic_Contrast);
        %%
    end
end
%%
[~,ind]=max(results_perc);
ind1 = logical(mod(ind,2));
ind2 = [~ind1(2:end) true];
figure(1);
plot(freqs,results_perc'); hold on; set(gca,'ColorOrderIndex',1);
plot(freqs(ind1),results_perc(1,ind1),'LineWidth',3); set(gca,'ColorOrderIndex',2);
plot(freqs(ind2),results_perc(2,ind2),'LineWidth',3); hold off;
xlim([1e2 1e4]);
grid minor
set(gca,'XScale','log');
title('Multizone vs Parametic - Percent Difference  (Third Octave)');
legend({'Orthogonal Basis Expansion Multizone';'Parametric Multizone'},'location','southwest');
ylabel('Percent Level Difference (%)');
xlabel('Frequency (Hz)');

[~,ind]=max(results_contrastdb);
ind1 = logical(mod(ind,2));
ind2 = [~ind1(2:end) true];
figure(2);
plot(freqs,results_contrastdb'); hold on; set(gca,'ColorOrderIndex',1);
plot(freqs(ind1),results_contrastdb(1,ind1),'LineWidth',3); set(gca,'ColorOrderIndex',2);
plot(freqs(ind2),results_contrastdb(2,ind2),'LineWidth',3); hold off;
xlim([1e2 1e4]);
grid minor
set(gca,'XScale','log');
title('Multizone vs Parametic - Acoustic Contrast  (Third Octave)');
legend({'Orthogonal Basis Expansion Multizone';'Parametric Multizone'},'location','southwest');
ylabel('Acoustic Contrast (dB)');
xlabel('Frequency (Hz)');

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script