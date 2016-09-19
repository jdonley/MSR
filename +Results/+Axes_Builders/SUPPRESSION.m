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
    0.2 0.2 1.0       ; ...       Known Signal
    1.0 0.0 0.0       ;]}; %      Predicted Signal

markers = {[ ...          Marker Shapes
    'o'       ; ...       Known Signal
    '>'       ;]}; %      Predicted Signal

trendAlpha = 0.5;

axes(axH);
[ax,h1,h2]=plotyy(1,1,1,1);delete(h1);delete(h2);

axs = ax;
set(axs,'Color','none'); %Transparent axes backgrounds
legendStrings = {};

%% Get results file paths
for rt = 1:numel(measures)
    Results_filepath = [ ...
        Results.getResultsPath( SYS  ), ...
        measures{rt}, '_Results.csv'];
    
    Res_Matrix = cell(2,1);
    Res_trend = cell(2,1);
    Res_area = cell(2,1);
    Res_CI = cell(2,1);
    
    %% Read results
    results_func = str2func(['Results.import_' results_types{rt} '_Reverb']);
    switch results_types{rt}          
        case 'Suppression'
            % Read results
            [NoiseLevel,Result_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] ...
                = results_func(Results_filepath);
            if NoiseLevel == -1
                Results.Axes_Builders.Helpers.setErrorAxis(ax,Result_Bright,Results_filepath)
                break
            end
            % Generate Plottable data matrices and vectors
            [Hrz_Vec, Res_Matrix{1}, Res_trend{1}, Res_area{1}, Res_CI{1}, CI_vec] = ...
                Results.generatePlotData( NoiseLevel, Result_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up, 'smoothingspline', [1.8 1.8]);            
            legendStrings = {legendStrings{:}, [measures{rt} ' BZ']};
            
        otherwise
            error(['Result type ''' results_types{rt} ''' not recognised']);
    end
    
    %%
    cols = colours{rt};
    mrks = markers{rt};
    
    domain = Hrz_Vec([1 end]);
    domain_lbl = 'Block Length ($\mathrm{ms}$)';
    
    switch results_types{rt}
        case 'Suppression'
            axCurr = ax(2);
            range = [-25 0];
            range_lbl = 'Suppression ($\mathrm{dB}$)';
    end
    
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
        
        erPl = errorbar(axCurr,Hrz_Vec,mean(Res_Matrix{pl}),Res_CI{pl}(:,1),Res_CI{pl}(:,2),...
            erPlDetails{:});
                
        hold off;
        
    end
    
end

grid(axs(2),'off'); % Turn off the second y axis grid because both y axes have the same grid

end



