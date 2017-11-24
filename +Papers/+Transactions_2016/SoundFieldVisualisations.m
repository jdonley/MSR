clc;
clear;
close all;
tic;

%%
% SYS = Current_Systems.loadCurrentSRsystem;
SYS = Current_Systems.IEEETransactions_System_G;

%%
setup = [SYS.Main_Setup(:);SYS.Masker_Setup(:)];
for s = 1:numel(setup)
    setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createSoundfield('DEBUG');
    setup(s) = setup(s).calc_Loudspeaker_Weights();
    setup(s) = setup(s).reproduceSoundfield('DEBUG');
end

%%
if exist('fH'), if isvalid(fH), close(fH); end; end


plotType = 'real';
% plotType = 'abs';


figNums = [101,102,103];
realistic = false;
details.PlotType = 'image'; % 'image' or 'surf'
details.DrawDetails = false;
details.zoneLineWid = 1.0;
details.arrowLineWid = 0.5;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 10;
details.PlotMic = false;
details.sphereRad = 2; % centimetres
details.NTicks = SYS.publication_info.axes_NumTicks;

for s = 1:numel(setup)
    %     reproZone = Orthogonal_Basis_Expansion.spatial_zone( ...
    %         setup(s).Multizone_Soundfield.k_global/2/pi*343, ...
    %         0, setup(s).Multizone_Soundfield.Radius, 'pw');
    %     reproZone.res = setup(s).res;
    %     reproZone = reproZone.createEmptySoundfield();
    %     reproRegionSamples = setup(s).getZoneSamples(reproZone);
    %     bMask = setup(s).Multizone_Soundfield.Bright_Zone.Soundfield_d_mean_mask;
    %     setup(s) = setup(s).save_Bright_Samples();
    %     maxBrightVal = mean(abs(setup(s).Bright_Samples(bMask(:))));
    %     reproRegionSamples = reproRegionSamples / maxBrightVal;
    %
    %     pk(s) = max(abs((reproRegionSamples(:))))*setup(s).res;
    pk(s) = mean(abs((setup(s).Bright_Samples(:))),'omitnan')*setup(s).res;
    %      pk(s) = max(abs((setup(s).Soundfield_reproduced(:))))*setup(s).res;
    
    F{s} = setup(s).Soundfield_reproduced*setup(s).res;
end


% clipFact = 2;
for s = 1:numel(setup)
    gainNorm = 1/max(pk(s)); % Normalise to the maximum of all bright peaks
    F{s} = F{s}*gainNorm; pk(s) = pk(s)*gainNorm;
    
    %     F{s}(abs(F{s})>clipFact*pk(s))=nan;
end


% close all;
fH = figure('Name',SYS.publication_info.FigureName);
fH.Color = 'w';
dimSz = SYS.publication_info.subPlotDims;
ha = tightPlots( ...
    dimSz(2), ...
    dimSz(1), ...
    SYS.publication_info.figure_width, ...
    SYS.publication_info.axis_aspect_ratio, ...
    SYS.publication_info.axes_gap, ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');
FontSize = 10;
FontName = 'Times';
FF = SYS.publication_info.LaTeX_FontFamily; % Set fonts for latex interpreter
FFnums = SYS.publication_info.LaTeX_NumbersFontFamily; % Set number fonts for latex interpreter
FS = SYS.publication_info.FontSize;
latexFontSettings = ['{\fontfamily{' FF '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];
latexNumFontSettings = ['{\fontfamily{' FFnums '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];


for s = 1:numel(setup)
    axes(ha(s));
    ax=gca;
    if contains(plotType,'abs')
        setup(s).plotSoundfield( mag2db(abs(F{s})), 'scientific_L12', realistic, details);
    elseif contains(plotType,'real')
        setup(s).plotSoundfield( real(F{s}), 'scientific_D1', realistic, details);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if contains(plotType,'abs')
        cm = cmap('L12','N',24);
        labcm = rgb2lab(cm);
        cmA = atan2(labcm(:,3),labcm(:,2))/pi*180;
        cmC = sum(labcm(:,2:3).^2,2).^.5;
        
        cmA(end-4:end,:) = cmA(end-4:end,:) + 45; % Red
        cmA(end-6:end-5,:) = cmA(end-6:end-5,:) + 270; % Green
        
        labcm(:,2:3) = cmC.*[cos(cmA/180*pi) sin(cmA/180*pi)];
        cm = lab2rgb(labcm);
        
        colormap(ax,cm);
    elseif contains(plotType,'real')
        cm = cmap('D1','N',24);
        colormap(ax,cm);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    txH = text(0,size(F{s},1),1e3,...
        ['(' char(64+s) '): ' ...
        SYS.publication_info.subPlotTitles{s}],...
        ... 'BackgroundColor',[1 1 1 0.7], ...
        'FontName',FontName,'FontSize',FontSize, ...
        'Interpreter',SYS.publication_info.Interpreter, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    if strcmpi(SYS.publication_info.Interpreter, 'latex')
        txH.String = [latexFontSettings txH.String '}'];
    end
    if s ~= 1
        ax.Title = [];
    else
        ax.Title.String = ' ';
    end
    
    ax.FontSize = SYS.publication_info.FontSize;
    ax.XLabel.FontSize = SYS.publication_info.FontSize;
    ax.YLabel.FontSize = SYS.publication_info.FontSize;
    
    ax.TickLabelInterpreter = SYS.publication_info.Interpreter;
    ax.XLabel.Interpreter = SYS.publication_info.Interpreter;
    ax.YLabel.Interpreter = SYS.publication_info.Interpreter;
    
    
    [c,r] = ind2sub(dimSz,s);
    if c ~= 1
        ax.YLabel = [];
        ax.YTickLabel = [];
    else
        if strcmpi(ax.TickLabelInterpreter, 'latex')
            for t = 1:numel(ax.YTickLabel)
                ax.YTickLabel(t) = {[latexNumFontSettings ax.YTickLabel{t} '}']};
            end
        end
    end
    if r ~= dimSz(2)
        ax.XLabel = [];
        ax.XTickLabel = [];
    else
        if strcmpi(ax.TickLabelInterpreter, 'latex')
            for t = 1:numel(ax.XTickLabel)
                ax.XTickLabel(t) = {[latexNumFontSettings ax.XTickLabel{t} '}']};
            end
        end
    end
    
    
    if contains(plotType,'real')
        clim_=[-1 1].*pk(s);
    elseif contains(plotType,'abs')
        clim_=[-18 6];
    end
    ax.CLim = clim_;
    colorbar off;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If bottom left axis then keep ylabel and shift to centre
for a = 1:dimSz(2)
    I = sub2ind(dimSz,1,a);
    axes(ha(I)); axY(a) = gca;
end

origUnits = get([axY axY(end).YLabel], 'Units');
set([axY axY(end).YLabel], 'Units', 'centimeters');

axY(end).YLabel.Position(2) = (sum(axY(1).Position([2 4])) - axY(end).Position(2))/2;
for a = 1:dimSz(2)-1
    axY(a).YLabel.String = '';
end

cellfun(@set, num2cell([axY axY(end).YLabel]'), ...
    repmat({'Units'},numel(origUnits),1), origUnits);
drawnow;

% If bottom right axis then keep xlabel and shift to centre
for a = 1:dimSz(1)
    I = sub2ind(dimSz,a,dimSz(2));
    axes(ha(I)); axX(dimSz(1)-a+1) = gca;
end
origUnits = get([axX axX(end).XLabel], 'Units');
set([axX axX(end).XLabel], 'Units', 'centimeters');

axX(end).XLabel.Position(1) = (sum(axX(1).Position([1 3])) - axX(end).Position(1))/2;
for a = 2:dimSz(1)
    axX(dimSz(1)-a+1).XLabel.String = '';
end

cellfun(@set, num2cell([axX axX(end).XLabel]'), ...
    repmat({'Units'},numel(origUnits),1), origUnits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if contains(plotType,'abs')
    cbTickSep = 3; %dB
elseif contains(plotType,'real')
    cbTickSep = 0.2; %Unitless
end
NcbTicks = numel(ax.CLim(1):cbTickSep:ax.CLim(end));
%
drawnow;
hCB = colorbar;
hCB.Location = 'manual';
hCB.Units = 'points';
hCB.TickDirection = SYS.publication_info.axes_tickdir;
hCB.Label.Interpreter = SYS.publication_info.Interpreter;
hCB.TickLabelInterpreter = SYS.publication_info.Interpreter;
hCB.Ticks = interp1(1:length(caxis),caxis,linspace(1,length(caxis),NcbTicks));
hCB.TickLabels = cellfun(@strrep, ...
    num2cell(num2str(linspace( ax.CLim(1), ax.CLim(2),NcbTicks)' ),2), ...
    repmat({' '},NcbTicks,1), repmat({''},NcbTicks,1),'un',0);
if strcmpi(hCB.TickLabelInterpreter, 'latex')
    hCB.TickLabels = cellfun(@(s) [latexNumFontSettings s '}'], hCB.TickLabels,'un',0);
end

if contains(plotType,'abs')
    hCB.Label.String = 'Magnitude (dB)';
elseif contains(plotType,'real')
    hCB.Label.String = 'Normalised Amplitude';
end

ax.Units = 'centimeters'; hCB.Units = 'centimeters';

hCB.Position(1) = ...
    ax.Position(1) + ...
    ax.Position(3) + ...
    SYS.publication_info.axes_gap(2);
hCB.Position(2) = ...
    ax.Position(2);
hCB.Position(3) = ...
    SYS.publication_info.colorbarWidth;
hCB.Position(4) = ...
    ax.Position(4) + ...
    SYS.publication_info.axes_gap(1) + ...
    ax.Position(4); % Colorbar height


tightfigadv;

savePath = [SYS.publication_info.DocumentPath, ...
    filesep, SYS.publication_info.FigureName, ...
    '.', SYS.publication_info.print_fmt];
% print(savePath,['-d' SYS.publication_info.print_fmt]);
export_fig(savePath,['-d' SYS.publication_info.print_fmt]);



%%
round([setup(1:2:end).k_global],SYS.publication_info.sigRounding,'s')
round([setup(1:2:end).k_global]/2/pi*SYS.signal_info.c,SYS.publication_info.sigRounding,'s')


%%
disp(['   Contrast: ' num2str(pow2db(setup(1).Acoustic_Contrast)) 'dB']);
disp(['        MSE: ' num2str(mag2db(setup(1).MSE_Bright)) 'dB']);
disp(['Attenuation: ' num2str(setup(1).Attenuation_dB(1)) 'dB (±' num2str(setup(1).Attenuation_dB(2)) 'dB)']);

%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script