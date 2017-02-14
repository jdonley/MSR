function [axs, legendStrings] = SUPPRESSION( SYS, axH )
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
% Copyright: Jacob Donley 2016
% Date: 04 September 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measures = SYS.analysis_info.Measures;
results_types = measures;

colours = {[ ...            R G B  values
    1.0 0.0 0.0       ; ...       Predicted Signal
    0.2 0.2 1.0       ;]}; %      Known Signal

markers = {[ ...          Marker Shapes
    'x'       ; ...       Predicted Signal
    '^'       ;]}; %      Known Signal

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
    results_func = str2func(['Results.import_' 'Suppression_Reverb']);
    
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
                Results.generatePlotData( Xvec, Vals(:,2), ConfInt_Low(:,2), ConfInt_Up(:,2), 'smoothingspline', [1.8 1.8]);
            Res_CI{1} = [ConfInt_Low(:,2), ConfInt_Up(:,2)];
            legendStrings = {legendStrings{:}, [measures{rt}]};
            
        case 'Actual Signal'
            % Generate Plottable data matrices and vectors
            [Hrz_Vec, Res_Matrix{2}, Res_trend{2}, Res_area{2}, Res_CI{2}, CI_vec] = ...
                Results.generatePlotData( Xvec, Vals(:,1), ConfInt_Low(:,1), ConfInt_Up(:,1), 'smoothingspline', [1.8 1.8]);
            Res_CI{2} = [ConfInt_Low(:,1), ConfInt_Up(:,1)];
            legendStrings = {legendStrings{:}, [measures{rt}]};
            
        otherwise
            error(['Result type ''' results_types{rt} ''' not recognised']);
    end
    
end
%%
cols = colours{:};
mrks = markers{:};

domain = Hrz_Vec([1 end]);
domain_lbl = 'Block Length ($\mathrm{ms}$)';

axCurr = ax(1);
range = [-20 0];
range_lbl = 'Suppression ($\mathrm{dB}$)';



Results.Axes_Builders.Helpers.setAxisParameters( SYS, axCurr, range, domain, range_lbl, domain_lbl);

for pl = size(cols,1):-1:1 % for each plottable set of data (plot in reverse order to keep axes order)
    
    PlDetails = { ...
        'Color',cols(pl,:), ...
        'LineWidth',SYS.publication_info.lineWid};
    erPlDetails = { ...
        mrks(pl), ...
        PlDetails{:}, ...
        'MarkerSize',4};
    
    axes(axCurr);
    hold on;
    
    trend_vec = linspace(Hrz_Vec(1),Hrz_Vec(end),numel(Hrz_Vec)*10);
    Pl = plot(axCurr,trend_vec,Res_trend{pl}(trend_vec),PlDetails{:});
    Pl.Color(4) = trendAlpha;
    
    erPl = errorbar(axCurr,Hrz_Vec, Res_Matrix{pl} , Res_CI{pl}(:,1), Res_CI{pl}(:,2),...
        erPlDetails{:});
    
    hold off;
    
end


% grid(axs(2),'off'); % Turn off the second y axis grid because both y axes have the same grid

end



