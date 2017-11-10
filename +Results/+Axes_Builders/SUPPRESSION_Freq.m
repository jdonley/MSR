function [axs, legendStrings] = SUPPRESSION_Freq( SYS, axH )
%SUPPRESSION Generates axes objects for suppression results
%
% Syntax:	[ axs, legendStrings ] = SUPPRESSION( SYS, axH )
%
% Inputs:
% 	SYS - Soundfield reproduction system specification object
%   axH - Axis handle to reuse
%
% Outputs:
% 	axs - Axis objects containing the STOI and PESQ styled plots.
%   legendEntries - The text for the legend of the plot.
%
%
% See also:

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016-2017
% Date: 4 September 2017
% Version: 0.2 (4 September 2017)
% Version: 0.1 (4 September 2016)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measures = SYS.analysis_info.Measures;
if isfield(SYS.analysis_info,'Result_Type') ...
        && ~isempty(SYS.analysis_info.Result_Type)
    results_types = SYS.analysis_info.Result_Type;
else
    results_types = measures;
end

colours = {[ ...            R G B  values
    1.0 0.0 0.0       ; ...       Predicted Signal
    0.2 0.2 1.0       ;]}; %      Known Signal

lineStyles = { ...          Marker Shapes
    '- '       ; ...       Predicted Signal
    '--'       ;}; %      Known Signal

markers = {[ ...          Marker Shapes
    '.'       ; ...       Predicted Signal
    '.'       ;]}; %      Known Signal

trendAlpha = 0.5;

axes(axH);
[ax,h1,h2]=plotyy(1,1,1,1);delete(h1);delete(h2);
ax(2).delete;ax(2)=[];

axs = ax;
set(axs,'Color','none'); %Transparent axes backgrounds
legendStrings = {};

%%
Res_Matrix = cell(numel(measures),1);
Res_trend = cell(numel(measures),1);
Res_area = cell(numel(measures),1);
Res_CI = cell(numel(measures),1);

for rt = 1:numel(measures)
    
    
    %% Read results
    results_func = str2func(['Results.import_' 'Suppression_vs_Freq']);
    
    % Read results
    [Xvec,Vals,ConfInt_Low,ConfInt_Up] ...
        = results_func( SYS );
    if Xvec == -1
        Results.Axes_Builders.Helpers.setErrorAxis(ax,Vals)
        break
    end
    
    switch results_types{rt}
        case 'Predicted Signal'
            % Generate Plottable data matrices and vectors
            [Hrz_Vec, Res_Matrix{1}, Res_trend{1}, Res_area{1}, Res_CI{1}, CI_vec] = ...
                Results.generatePlotData( Xvec, Vals(:,1), ConfInt_Low(:,1), ConfInt_Up(:,1), 'smoothingspline', [1.8 1.8]);
            Res_CI{1} = [ConfInt_Low(:,1), ConfInt_Up(:,1)];
            legendStrings = {legendStrings{:}, [results_types{rt}]};
            
        case 'Actual Signal'
            % Generate Plottable data matrices and vectors
            [Hrz_Vec, Res_Matrix{2}, Res_trend{2}, Res_area{2}, Res_CI{2}, CI_vec] = ...
                Results.generatePlotData( Xvec, Vals(:,2), ConfInt_Low(:,2), ConfInt_Up(:,2), 'smoothingspline', [1.8 1.8]);
            Res_CI{2} = [ConfInt_Low(:,2), ConfInt_Up(:,2)];
            legendStrings = {legendStrings{:}, [results_types{rt}]};
            
        otherwise
            error(['Result type ''' results_types{rt} ''' not recognised']);
    end
    
end
%%
cols = colours{:};
lines= lineStyles(:);
mrks = markers{:};

domain = round(Hrz_Vec([1 end]),SYS.publication_info.sigRounding,'significant');
domain_lbl = 'Frequency ($\mathrm{kHz}$)';

axCurr = ax(1);
range = [-15 5];
range_lbl = 'Suppression ($\mathrm{dB}$)';



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

Falias = 2;

plA = plot([1 1]*Falias,range*2,'--k');
plA.Tag = 'noLegend';
    
ha=area(repmat([Falias,ceil(Hrz_Vec(end))],2,1).', ...
    repmat([range(1)*2 diff(range*2)],2,1),'LineStyle','none','showbaseline','off');
ha(1).FaceColor=[1 1 1];drawnow;
ha(2).FaceColor=[0 0 0];drawnow;
if V >= 8.6 % If MATLAB version is greater than 8.6 (R2015b)
    ha(1).FaceAlpha=0.0;
    ha(2).FaceAlpha=0.2;
else
    ha(1).Face.ColorType = 'truecoloralpha';
    ha(1).Face.ColorData(4) = 255 * 0.0;
    ha(2).Face.ColorType = 'truecoloralpha';
    ha(2).Face.ColorData(4) = 255 * 0.2;
end


plB = plot([1 1]*domain(1),range*2,'-k');
plB.Tag = 'noLegend';

ha=area(repmat([0.01,domain(1)],2,1).', ...
    repmat([range(1)*2 diff(range*2)],2,1),'LineStyle','none','showbaseline','off');
ha(1).FaceColor=[1 1 1];drawnow;
ha(2).FaceColor=[0 0 0];drawnow;

hold off;

end



