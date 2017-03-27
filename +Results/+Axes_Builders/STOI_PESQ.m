function [axs, legendStrings] = STOI_PESQ( SYS, axH )
%STOI_PESQ Generates axes objects with STOI and PESQ results
%
% Syntax:	[ ax ] = STOI_PESQ( SYS )
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
% Date: 09 June 2016
% Revision: 0.2 (23 March 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measures = {'STOI';'PESQ'};
results_types = {'SpeechIntelligibility';'Quality'};

colours = {[ ...            R G B  values
    0.2 0.2 1.0       ; ...       Bright Intelligibility
    1.0 0.0 0.0       ;];[ ...    Quiet  Intelligibility
    0.6 0.0 0.6       ];}; %      Bright Quality

markers = {[ ...          Marker Shapes
    'o'       ; ...       Bright Intelligibility
    '>'       ;];[ ...    Quiet  Intelligibility
    'd'       ];}; %      Bright Quality

lineStys = { ...         Line Styles
    '-'      ; ...       Line 1
    '--'     ; ...       Line 2
    ':'      }; %        Line 3

trendAlpha = 0.5;

axes(axH);
[ax,h1,h2]=plotyy(1,1,1,1);delete(h1);delete(h2);

axs = ax;
set(axs,'Color','none'); %Transparent axes backgrounds
legendStrings = {};

mergeLines = SYS.publication_info.MergeLines;


%% Repeat for each line on axis
SYS_AllLines = SYS;
Nlines = numel(SYS.Main_Setup);

if mergeLines
   lineStys{Nlines} = '-';
end

for li = 1:Nlines
    SYS.Main_Setup = SYS_AllLines.Main_Setup(li);
    SYS.Masker_Setup = SYS_AllLines.Masker_Setup(li);
    
    %% Get results file paths
    for rt = 1:numel(measures)
        Results_filepath = [ ...
            Results.getResultsPath( SYS  ), ...
            measures{rt}, '_Results.csv'];
        
        Res_Matrix  = cell(2,1);
        Res_trend   = cell(2,1);
        Res_area    = cell(2,1);
        Res_CI      = cell(2,1);        
        
        %% Read results
        results_func = str2func(['Results.import_' results_types{rt} '_Reverb']);
        switch results_types{rt}
            case 'SpeechIntelligibility'
                % Read results
                [NoiseLevel,Result_Bright,Result_Quiet,ConfInt_Bright_Low,ConfInt_Bright_Up,ConfInt_Quiet_Low,ConfInt_Quiet_Up] ...
                    = results_func(Results_filepath);
                if NoiseLevel == -1
                    Results.Axes_Builders.Helpers.setErrorAxis(ax,Result_Bright,Results_filepath)
                    break
                end
                % Remove values outside plot area according to system settings
                rmI = NoiseLevel<min(SYS.signal_info.L_noise_mask) | NoiseLevel>max(SYS.signal_info.L_noise_mask);
                NoiseLevel(rmI)=[];Result_Bright(rmI)=[];Result_Quiet(rmI)=[];
                ConfInt_Bright_Low(rmI)=[];ConfInt_Bright_Up(rmI)=[];ConfInt_Quiet_Low(rmI)=[];ConfInt_Quiet_Up(rmI)=[];
                
                % Generate Plottable data matrices and vectors
                [Hrz_Vec, Res_Matrix{1}, Res_trend{1}, Res_area{1}, Res_CI{1}, CI_vec] ...
                    = Results.generatePlotData( NoiseLevel, Result_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up);
                if li == 1
                    legendStrings = {legendStrings{:}, [measures{rt} ' BZ']};
                end
                [~, Res_Matrix{2}, Res_trend{2}, Res_area{2}, Res_CI{2}, ~] ...
                    = Results.generatePlotData( NoiseLevel, Result_Quiet, ConfInt_Quiet_Low, ConfInt_Quiet_Up);
                if li == 1
                    legendStrings = {legendStrings{:}, [measures{rt} ' QZ']};
                end
                
            case 'Quality'
                % Read results
                [NoiseLevel,Result_Bright,ConfInt_Bright_Low,ConfInt_Bright_Up] ...
                    = results_func(Results_filepath);
                if NoiseLevel == -1
                    Results.Axes_Builders.Helpers.setErrorAxis(ax,Result_Bright,Results_filepath)
                    break
                end
                % Remove values outside plot area according to system settings
                rmI = NoiseLevel<min(SYS.signal_info.L_noise_mask) | NoiseLevel>max(SYS.signal_info.L_noise_mask);
                NoiseLevel(rmI)=[];Result_Bright(rmI)=[];
                ConfInt_Bright_Low(rmI)=[];ConfInt_Bright_Up(rmI)=[];
                
                % Generate Plottable data matrices and vectors
                [Hrz_Vec, Res_Matrix{1}, Res_trend{1}, Res_area{1}, Res_CI{1}, CI_vec] = ...
                    Results.generatePlotData( NoiseLevel, Result_Bright, ConfInt_Bright_Low, ConfInt_Bright_Up, 'smoothingspline', [1.8 1.8]);
                if li == 1
                    legendStrings = {legendStrings{:}, [measures{rt} ' BZ']};
                end
                
            otherwise
                error(['Result type ''' results_types{rt} ''' not recognised']);
        end
        
        %%
        trend_vec = linspace(Hrz_Vec(1),Hrz_Vec(end),numel(Hrz_Vec)*10);
        cfArgs={'UniformOutput', false};
        if mergeLines && li>1
            
            Res_trend_{rt}  = cellfun( @(v1,v2) v1(trend_vec)+v2, ...
                Res_trend(~cellfun('isempty',Res_trend)), ...
                Res_trend_{rt}(~cellfun('isempty',Res_trend_{rt})), 'un',0);
            
            Res_Matrix_{rt} = cellfun( @(v1,v2) v1 + v2 , ...
                Res_Matrix(~cellfun('isempty',Res_Matrix)), ...
                Res_Matrix_{rt}(~cellfun('isempty',Res_Matrix_{rt})), 'un',0);
            
            Res_CI_{rt}     = cellfun( @(v1,v2) v1 + v2 , ...
                Res_CI(~cellfun('isempty',Res_CI)), ...
                Res_CI_{rt}(~cellfun('isempty',Res_CI_{rt})), 'un',0);
            
        else            
            Res_trend_{rt}  = cellfun(@(v1) v1(trend_vec), ...
                Res_trend(~cellfun('isempty',Res_trend)), 'un',0);
            Res_Matrix_{rt} = Res_Matrix;
            Res_CI_{rt}     = Res_CI;
        end
        %%
        if (mergeLines && li==Nlines) || ~mergeLines
            if mergeLines
                Res_trend_{rt}  = cellfun(@(v1) v1/Nlines, ...
                    Res_trend_{rt}, 'un',0);
                Res_Matrix_{rt}  = cellfun(@(v1) v1/Nlines, ...
                    Res_Matrix_{rt}, 'un',0);
                Res_CI_{rt}  = cellfun(@(v1) v1/Nlines, ...
                    Res_CI_{rt}, 'un',0);
            end
            
            
            cols = colours{rt};
            mrks = markers{rt};
            
            domain = Hrz_Vec([1 end]);
            domain_lbl = 'Noise Mask Level ($G$) ($\mathrm{dB}$)';
            
            switch results_types{rt}
                case 'SpeechIntelligibility'
                    axCurr = ax(1);
                    range = [0 100];
                    range_lbl = 'STOI (\%WC)';
                case 'Quality'
                    axCurr = ax(2);
                    range = [1 4.56];
                    range_lbl = 'PESQ (MOS)';
            end
            
            Results.Axes_Builders.Helpers.setAxisParameters( SYS, axCurr, range, domain, range_lbl, domain_lbl);
            
            for pl = size(cols,1):-1:1 % for each plottable set of data (plot in reverse order to keep axes order)
                
                PlDetails = { ...
                    'Color',cols(pl,:), ...
                    'LineWidth',SYS.publication_info.lineWid};
                
                if Nlines > 1
                    PlOnlyDetails = {'LineStyle', lineStys{li}};
                else
                    PlOnlyDetails={};
                end
                
                erPlDetails = { ...
                    mrks(pl), ...
                    PlDetails{:}, ...
                    'MarkerSize',SYS.publication_info.markerSize};
                
                v=ver('matlab');
                if str2num(v.Version)>=9.2
                    erPlDetails = {erPlDetails{:}, 'CapSize',SYS.publication_info.capSize};
                end
                
                axes(axCurr);
                hold on;
                
                Pl = plot(axCurr,trend_vec,Res_trend_{rt}{pl},PlDetails{:},PlOnlyDetails{:});
                Pl.Color(4) = trendAlpha;
                
                erPl = errorbar(axCurr,Hrz_Vec,mean(Res_Matrix_{rt}{pl}),Res_CI_{rt}{pl}(:,1),Res_CI_{rt}{pl}(:,2),...
                    erPlDetails{:});
                
                hold off;
                
                if li ~= 1 && ~SYS.publication_info.MergeLines
                    Pl.Tag = 'noLegend';
                    erPl.Tag = 'noLegend';
                end
                
            end
            
        end
        
    end
    
end

grid(axs(2),'off'); % Turn off the second y axis grid because both y axes have the same grid
axs(2).Box = 'off';

end



