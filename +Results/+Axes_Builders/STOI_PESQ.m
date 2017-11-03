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
    1.0 0.0 0.0       ; ...       Quiet  Intelligibility
    0.6 0.0 0.6       ;];[ ...    Speech Intelligibility Contrast
    0.0 0.6 0.0       ];}; %      Bright Quality

markers = {[ ...          Marker Shapes
    'o'       ; ...       Bright Intelligibility
    '>'       ; ...       Quiet  Intelligibility
    's'       ;];[ ...    Speech Intelligibility Contrast
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

if isfield(SYS.publication_info,'MergeLines')
    mergeLines = SYS.publication_info.MergeLines;
else
    mergeLines = false;
end

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
                %%% Include Speech Intelligibility Contrast (SIC) in the plot
                if li == 1
                    legendStrings = {legendStrings{:}, ['SIC']};
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
            
            Res_Matrix_{rt} = cellfun( @(v1,v2) [v1; v2] , ...
                Res_Matrix(~cellfun('isempty',Res_Matrix)), ...
                Res_Matrix_{rt}(~cellfun('isempty',Res_Matrix_{rt})), 'un',0);
            
%             Res_CI_{rt}     = cellfun( @(v1,v2) v1 + v2 , ...
%                 Res_CI(~cellfun('isempty',Res_CI)), ...
%                 Res_CI_{rt}(~cellfun('isempty',Res_CI_{rt})), 'un',0);
            
        else            
            Res_trend_{rt}  = cellfun(@(v1) v1(trend_vec), ...
                Res_trend(~cellfun('isempty',Res_trend)), 'un',0);
            Res_Matrix_{rt} = Res_Matrix;
            Res_CI_{rt}     = Res_CI;
        end
        
        if strcmpi(results_types{rt},'SpeechIntelligibility') && li==Nlines
            %%% Include Speech Intelligibility Contrast (SIC) in the plot
            Res_Matrix_{rt}{3} = Res_Matrix_{rt}{1} - Res_Matrix_{rt}{2};
            Res_trend_{rt}{3}  = Res_trend_{rt}{1}  - Res_trend_{rt}{2};
        end
        %%
        if (mergeLines && li==Nlines) || ~mergeLines
            if mergeLines
                Res_trend_{rt}  = cellfun(@(v1) v1/Nlines, ...
                    Res_trend_{rt}, 'un',0);
%                 Res_Matrix_{rt}  = cellfun(@(v1) v1/Nlines, ...
%                     Res_Matrix_{rt}, 'un',0);
%                 Res_CI_{rt}  = cellfun(@(v1) v1/Nlines, ...
%                     Res_CI_{rt}, 'un',0);
                Res_CI_{rt}  = cellfun(@(v1) Tools.confidence_intervals(v1,95), ... % 95 percent confidence interval
                    Res_Matrix_{rt}, 'un',0);
            end
            

            
            cols = colours{rt};
            mrks = markers{rt};
            
            domain = Hrz_Vec([1 end]);
            if isfield(SYS.publication_info,'xlabel')
                domain_lbl = SYS.publication_info.xlabel;
            else
                domain_lbl = 'Noise Mask Level ($G$) ($\mathrm{dB}$)';
            end
            
            switch results_types{rt}
                case 'SpeechIntelligibility'
                    axCurr = ax(1);
                    range = [0 100];
                    if isfield(SYS.publication_info,'ylabel')
                        range_lbl = SYS.publication_info.ylabel;
                    else
                        range_lbl = 'STOI (\%WC) or SIC (\%)';
                    end                    
                    % Confidential Speech Privacy Area Shading (WC < 25%)
                    % Determined from: ASTM E1130 and "ASTM METRICS FOR RATING SPEECH PRIVACY OF CLOSED ROOMS AND OPEN PLAN SPACES"
                    STOI_Conflvl = 25;
                    axes(axCurr); hold on;
                    %%%
                    tx = text(mean(domain),(STOI_Conflvl-0)/2,upper('confidential'));
                    tx.BackgroundColor = 'none';
                    tx.Color = (1 - (1 - colours{rt}(2,:))*0.2);
                    tx.HorizontalAlignment = 'center';
                    tx.FontWeight = 'bold';
                    %%%
                    arS = area(axCurr, ...
                        ([1;1]*domain*2)', ...
                        [-10*[1 1]; 10 + STOI_Conflvl*[1 1]]'); hold off;
                    set(arS,'FaceColor', colours{rt}(2,:),'FaceAlpha', 0.05,'EdgeColor', 'none','BaseValue', -10);
                    arS(1).BaseLine.Color = 'none';arS(1).FaceAlpha = 0;
                    %%%
                    % Good Speech Quality (PESQ MOS-LQO = 4 (3.5 to 4.5))
                    PESQ_Goodlvl = 4.0;
%                     axes(axCurr);
                    hold on;
                    %%%
                    tx = text(mean(domain),((PESQ_Goodlvl + 0.5)-1)/3.56*100,upper('Good Quality'));
                    tx.BackgroundColor = 'none';
                    tx.Color = (1 - (1 - colours{contains(lower(results_types),'quality')}(1,:))*0.2);
                    tx.HorizontalAlignment = 'center';
                    tx.FontWeight = 'bold';
                    %%%
                    arS = area(axCurr, ...
                        ([1;1]*domain*2)', ...
                        ([PESQ_Goodlvl*[1 1] - 1; 5.5 5.5]')/3.56*100); hold off;
                    set(arS,'FaceColor', colours{contains(lower(results_types),'quality')}(1,:),...
                        'FaceAlpha', 0.05,'EdgeColor', 'none','BaseValue', ((PESQ_Goodlvl)-1)/3.56*100);
                    arS(1).BaseLine.Color = 'none';arS(1).FaceAlpha = 0;
                    %%%
                    
                case 'Quality'
                    axCurr = ax(2);
                    range = [1 4.56];
                    if isfield(SYS.publication_info,'ylabel2')
                        range_lbl = SYS.publication_info.ylabel2;
                    else
                        range_lbl = 'PESQ (MOS-LQO) (WB)';
                    end

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
                
                if li ~= 1 && ~mergeLines
                    Pl.Tag = 'noLegend';
                    erPl.Tag = 'noLegend';
                end
                
            end
            
        end
        
    end
    
end

grid(axs(2),'off'); % Turn off the second y axis grid because both y axes have the same grid
axs(2).Box = 'off';

% Some testing analysis if needed
if rt==2
    G = trend_vec;
    STOIB=(Res_trend_{1}{1})/100;
    STOIQ=(Res_trend_{1}{2})/100;
    PESQB=(Res_trend_{2}{1}-1)/3.56;
    %                 sum(PESQB)-sum(STOIQ);
    SIC = STOIB-STOIQ;
    fHdetail = figure(12321);
    fHdetail.Position(1:2) = [100 100];
    txtOffs = [-0.5 +0.5 -0.5];
    lambdas = [0.33 1 3];
    for l = 1:numel(lambdas)
        fHdetail = figure(12321);
        lambda = lambdas(l);
        optCurve = SIC + lambda*PESQB;
        plot(G,optCurve); hold on;
        Iopt = optCurve==max(optCurve);
        Gopt = G(Iopt);
        plot(Gopt,max(optCurve),'or');
        text(Gopt,max(optCurve)+txtOffs(l),...
            {['\lambda = ' num2str(lambdas(l),3)]; ...
             ['G = ' num2str(Gopt,3)]; ...
             ['SIC_{STOI} = ' num2str(SIC(Iopt)*100,3) '%']; ...
             ['B_{PESQ} = ' num2str(PESQB(Iopt)*3.56+1,3) 'MOS']},'ho','c');
         
         if lambdas(l) == 1
             axes(axs(1)); hold on;
             plot(axs(1),Gopt*[1 1],axs(1).YLim,':k');
%              plot(axs(1),G,SIC*100,'-','color',[0,0,0,0.5]);
             hold off;
         end
    end
    hold off;
end
0; % Breakpoint here to stop at each plot
end

% function saveSelectResultsAsLATEXmacros
%    
% %
% Results = [Eps_min_dB, ...
%            Eps_minB_dB, ...
%            mag2db(mean(db2mag([Eps_min_dB;Eps_minB_dB])))];
% newcom = '\\newcommand';
% mac_pre = '\\';
% mac_post = 'res';
% mac_names = {'LTASS'; ...
%              'QuietOne'; ...
%              'LowPass'; ...
%              'White'; ...
%              'Pink'; ...
%              'IBLambdaZero'; ...
%              'IBLambdaHalf'; ...
%              'IBLambdaOne'; };
% mac_rownames = {''; ...
%                 'B'; ...
%                 'AVG';};
%             
% cols = numel(mac_names);
% rows = numel(mac_rownames);
% mac_post = '';
% mac_names = num2cell(char((1:cols).' +64));
% mac_rownames = num2cell(char((1:rows).' +64));           
% M = cellfun(@(a,b) [a,b],...
%     repmat({ mac_pre   },cols,1),...
%              mac_names ,...
%     'UniformOutput',false);
% macros = cellfun(@(a,b,c) [a,b,c],...
%     repmat( M ,rows,1),...
%     reshape(repmat( mac_rownames.' ,cols,1),[],1),...
%     repmat({ mac_post  },rows*cols,1),...
%     'UniformOutput',false);
%       
% FC1 = cellfun(@(m) [newcom '{' m '} {{' ],macros,'UniformOutput',false);
% FC2 = [num2str(round(Results.',3,'significant')) repmat('}}\n',numel(macros),1)];
% 
% filecontent = cellfun(@(a,b) [a,b],FC1,mat2cell(FC2,ones(size(FC2,1),1),size(FC2,2)),'UniformOutput',false);
% 
% fid = fopen([SYS.publication_info.DocumentPath filesep SYS.publication_info.LatexMacrosFile],'w');
% fprintf(fid,[filecontent{:}]);
% 
% fclose(fid);
% 
% Tools.MiKTeX_FNDB_Refresh;
% end


