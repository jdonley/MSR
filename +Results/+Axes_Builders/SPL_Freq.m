function [axs, legendStrings] = SPL_Freq( SYS, axH )
%SPL_Freq Generates axes objects for suppression results
%
% Syntax:	[ axs, legendStrings ] = SPL_Freq( SYS, axH )
%
% Inputs:
% 	SYS - Soundfield reproduction system specification object
%   axH - Axis handle to reuse
%
% Outputs:
% 	axs - Axis objects containing the SPL vs Frequency styled plots.
%   legendEntries - The text for the legend of the plot.
%
%
% See also:

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 18 January 2017
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigNames = SYS.analysis_info.SignalNames;
recTypes = SYS.signal_info.recording_type;

% Hyphenate "Realworld" word if it is not already
recTypes = regexprep(lower(recTypes),'realworld','real-world');

% Force capitalisation
recTypes = regexprep(lower(recTypes),'(\<[a-z])','${upper($1)}');


colours = {[ ...            R G B  values
    0.2 0.2 1.0       ; ...       Simulated Bright Zone
    1.0 0.0 0.0       ; ...       Simulated Quiet Zone
    0.4 0.0 0.8       ; ...       Realworld Bright Zone
    0.6 0.4 0.0       ;]}; %      Realworld Quiet Zone

lineStyles = { ...         Lines
    '- '       ; ...       Simulated Bright Zone
    ': '       ; ...       Simulated Quiet Zone
    '-.'       ; ...       Realworld Bright Zone
    '--'       ;}; %       Realworld Quiet Zone

markers = {[ ...          Marker Shapes
    '.'       ; ...       Simulated Bright Zone
    '.'       ; ...       Simulated Quiet Zone
    '.'       ; ...       Realworld Bright Zone
    '.'       ;]}; %      Realworld Quiet Zone

trendAlpha = 0.33;

axes(axH);
[ax,h1,h2]=plotyy(1,1,1,1);delete(h1);delete(h2);
ax(2).Visible = 'off';
ax(2).Title.Visible = 'on';
% ax(2).delete;ax(2)=[]; % Uncomment and delete these to remove the lettered labels from the subplots


axs = ax;
set(axs,'Color','none'); %Transparent axes backgrounds
legendStrings = {};

%%
RT = numel(recTypes);
SN = numel(sigNames);
Res_Matrix  = cell(RT,SN);
Res_CI      = cell(RT,SN);
    
    
    %% Read results
    results_func = str2func(['Results.import_' 'SoundPressureLevels_vs_Freq']);
    
    % Read results
    [Xvec,Vals,ConfInt_Low,ConfInt_Up] ...
        = results_func( SYS );
    if Xvec == -1
        Results.Axes_Builders.Helpers.setErrorAxis(ax,Vals{3})
        return
    end
    
    for rt = 1:RT
        for sn = 1:SN
            if any(rt==1:size(Vals,2)) && any(sn==1:size(Vals,3))
                % Generate Plottable data matrices and vectors
                [Hrz_Vec, Res_Matrix{rt,sn}] = ...
                    Results.generatePlotData( Xvec, Vals(:,rt,sn), ConfInt_Low(:,rt,sn), ConfInt_Up(:,rt,sn), 'smoothingspline', [1.8 1.8]);
                Res_CI{rt,sn} = [ConfInt_Low(:,rt,sn), ConfInt_Up(:,rt,sn)];
                legendStrings = {legendStrings{:}, [recTypes{rt} ' ' sigNames{sn}]};
            else
                Res_Matrix{rt,sn} = NaN(1,numel(Hrz_Vec));
                Res_CI{rt,sn}     = NaN(numel(Hrz_Vec),2);
                legendStrings = {legendStrings{:}, ['Not Found']};
            end
        end
    end
    Res_Matrix = Res_Matrix.';
    Res_CI     = Res_CI.';
    
%%
cols = colours{:};
lines= lineStyles(:);
mrks = markers{:};

domain = round(Hrz_Vec([1 end]),SYS.publication_info.sigRounding,'significant');
domain_lbl = 'Frequency ($\mathrm{kHz}$)';

axCurr = ax(1);
range = [-40 0];
if isfield(SYS.publication_info,'YLim_override')
    range = SYS.publication_info.YLim_override;
end
range_lbl = 'Magnitude ($\mathrm{dB}$)';



Results.Axes_Builders.Helpers.setAxisParameters( SYS, axCurr, range, domain, range_lbl, domain_lbl);

for pl = size(cols,1):-1:1 % for each plottable set of data (plot in reverse order to keep axes order)

    PlDetailsTrend = { ...
        'Color',cols(pl,:), ...
        'LineStyle', lines{pl}, ...
        'LineWidth',SYS.publication_info.lineWid };
    PlDetailsShade = { ...
        'LineStyle', 'none', ...
        'showbaseline','off'};    
    PlDetailsMark = {...
        'Color',cols(pl,:), ...
        'LineStyle', lines{pl}, ...
        'Marker', mrks(pl), ...
        'MarkerSize',5, ...
        'LineWidth',SYS.publication_info.lineWid };
    
    axes(axCurr);
    hold on;
    
    ha=area([Hrz_Vec; Hrz_Vec].', [Res_Matrix{pl}.'+Res_CI{pl}(:,2), -sum(Res_CI{pl},2)], ...
        PlDetailsShade{:} );
    ha(1).FaceColor=cols(pl,:);drawnow;
    ha(2).FaceColor=cols(pl,:);drawnow;
    v=version;i=find(v=='.');V=str2num(v(1:(i(2)-1)));
    if V >= 8.6 % If MATLAB version is greater than 8.6 (R2015b)
        ha(1).FaceAlpha=0.0;
        ha(2).FaceAlpha=trendAlpha;
    else
       ha(1).Face.ColorType = 'truecoloralpha';
       ha(1).Face.ColorData(4) = 255 * 0.0;
       ha(2).Face.ColorType = 'truecoloralpha';
       ha(2).Face.ColorData(4) = 255 * trendAlpha;
    end        
    
%     trend_vec = linspace(Hrz_Vec(1),Hrz_Vec(end),numel(Hrz_Vec)*10);
%     Pl = plot(axCurr,trend_vec,Res_trend{pl}(trend_vec),PlDetailsTrend{:});
    
    PlM = plot(axCurr,Hrz_Vec,Res_Matrix{pl},PlDetailsMark{:});
    

    

    

    
end

Falias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c / 1e3;

shadeRange = range+[-1 1]*diff(range);
plA = plot([1 1]*Falias,shadeRange,'--k');
plA.Tag = 'noLegend';
    
% ha=area(repmat([Falias,ceil(Hrz_Vec(end))],2,1).', ...
%     repmat([shadeRange(1) diff(shadeRange)],2,1),'LineStyle','none','showbaseline','off');
% ha(1).FaceColor=[1 1 1];drawnow;
% ha(2).FaceColor=[0 0 0];drawnow;
% if V >= 8.6 % If MATLAB version is greater than 8.6 (R2015b)
%     ha(1).FaceAlpha=0.0;
%     ha(2).FaceAlpha=0.2;
% else
%     ha(1).Face.ColorType = 'truecoloralpha';
%     ha(1).Face.ColorData(4) = 255 * 0.0;
%     ha(2).Face.ColorType = 'truecoloralpha';
%     ha(2).Face.ColorData(4) = 255 * 0.2;
% end


% plB = plot([1 1]*domain(1),range*2,'-k');
% plB.Tag = 'noLegend';

% ha=area(repmat([0.01,domain(1)],2,1).', ...
%     repmat([range(1)*2 diff(range*2)],2,1),'LineStyle','none','showbaseline','off');
% ha(1).FaceColor=[1 1 1];drawnow;
% ha(2).FaceColor=[0 0 0];drawnow;

hold off;

end



