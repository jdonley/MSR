%clc;
clear;
close all;
tic;

%%
DocumentPath = 'tex\latex\MSR_Setup_Diagram';
print_fmt = 'epsc'; %figure image file format
fileext = 'eps';
print_res = 600; %dpi
plot_width = 88.9/10;% + 6.35/10 + 88.9/10; %IEEE full text width
aspect_ratio = 4.2/3;
FontSize = 9;
Font = 'Times';

%%
room_width = 3.200; % metres
rad = ( room_width - 2*0.115 - 2*0.085 ) / 2; % metres


speech_layout = {'brightzone_pos_angle',        90, ...
    'quietzone_pos_angle',         -90, ...
    'brightzone_source_angle',     0, ...
    'brightzone_source_type',      'pw'};

masker_layout = {'brightzone_pos_angle',        -90, ...
    'quietzone_pos_angle',         90, ...
    'brightzone_source_angle',     180, ...
    'brightzone_source_type',      'ps'};

array_type = 'circle';
spkr_radius = 1.3;
spkr=1;

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
if spkr ==1
    %regular
    N_spkrs = 24;
    array = { ...
        'angleto_firstloudspeaker',      90, ...
        'angleof_loudspeakerarc',        180 * N_spkrs/(N_spkrs-1) , ...
        'numberof_loudspeakers',         N_spkrs, ...
        'loudspeaker_model',             'Genelec 8010A', ...
        'loudspeaker_radius',            spkr_radius, ...
        'loudspeaker_spacing',           spkr_spacing, ...
        'speaker_array_type',            array_type};
elseif spkr == 2
    %parametric
    array = { ...
        'angleto_firstloudspeaker',     th/pi*180, ...
        'loudspeaker_radius',           x_, ...
        'numberof_loudspeakers',        1, ...
        'loudspeaker_model',            'Parametric', ...
        'loudspeaker_spacing',          0.01, ...
        'speaker_array_type',           'line'};
end

setup = Speaker_Setup.createSetup({...
    'frequency',                    1000, ...
    speech_layout{:}, ...
    array{:}, ...
    'resolution',                   100, ... % Minimum resolution of approx 50 for 8kHz signal to satisfy nyquist theorem. We choose 100 for good measure.
    'reproduction_radius',          1.0, ...
    'bright_weight',                1.0, ...
    'quiet_weight',                 0, ...
    'unattended_weight',            0.05, ...
    'brightzone_radius',            0.3, ...
    'brightzone_source_dist',       x_, ...
    'brightzone_pos_distance',      0.6, ...
    'quietzone_radius',             0.3, ...
    'quietzone_pos_distance',       0.6, ...
    'maximum_frequency',            8000, ...
    'angleof_loudspeakerarrcentre', 180, ...
    'loudspeaker_object',           Parametric_Synthesis.parametric_soundfield});

%%
setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');

setup = setup.calc_Loudspeaker_Weights();
setup = setup.reproduceSoundfield('DEBUG');

%%
figure(123);

set(gcf, 'Units','centimeters', 'Color','w');
fig_pos = get(gcf,'Position');
set(gcf, 'PaperUnits','centimeters', ...
    'Position', [fig_pos(1) fig_pos(2) plot_width plot_width/aspect_ratio], ...
    'PaperSize', [plot_width plot_width/aspect_ratio]);

hold off;
realistic = true;
details.DrawDetails = true;
details.zoneLineWid = 1.5;
details.arrowLineWid = 0.4;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 30;

pk = max(abs(real(setup.Bright_Samples(:))));
Z = setup.Soundfield_reproduced;

setup.plotSoundfield( (Z), 'default', realistic, details);
plH = gca;

if spkr ==1
    R = [0 2].*spkr_radius*100 ; xlim(R);ylim(R);
elseif spkr == 2
    R = [0 2].*spkr_radius*100 - (spkr_radius-x_)*100; xlim(R);ylim(R);
end
%caxis([-pk, pk] );


%title('Small Zone Weight');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustments for publication
lims = [xlim(plH) ylim(plH) zlim(plH)];
axis(plH, lims+[-1 1 -1 1 -0.2 2]*15);

waveFieldAmplification = 0.01;
waveFieldShift = -1;
hTmp=findobj(plH.Children,'Type','Surface');
hTmp(end).ZData = hTmp(end).ZData.*waveFieldAmplification+waveFieldShift;
plH.CLim = plH.CLim.*waveFieldAmplification+waveFieldShift;

colorbar off
grid on
plH.Title=[];
plH.XLabel=[];
plH.YLabel=[];
plH.ZLabel=[];
plH.XTickLabel=[];
plH.YTickLabel=[];
plH.ZTickLabel=[];
plH.TickLength=[0 0];


camproj('perspective');
view([90 35]);
camzoom(2.15);


%% Save Figure
if ~exist(DocumentPath,'dir'); mkdir(DocumentPath); end
%export_fig([DocumentPath '\MSR_Layout.' fileext], ['-r' num2str(print_res)]);
%print([DocumentPath '\MSR_Layout_' array_type '.' fileext], ['-d' print_fmt], ['-r' num2str(print_res)]);

%close all;
% Update latex File Name DataBase
%Tools.MiKTeX_FNDB_Refresh;

%%
fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script