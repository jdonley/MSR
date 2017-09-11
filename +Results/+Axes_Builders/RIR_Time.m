function [axs, legendStrings] = RIR_Time( SYS, axH )
%SPL_Freq Generates axes objects for RIR results
%
% Syntax:	[ axs, legendStrings ] = RIR_Time( SYS, axH )
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
% Date: 10 September 2017
% Version: 0.1 (10 September 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigNames = SYS.analysis_info.SignalNames;
% recTypes = SYS.signal_info.recording_type;

% Hyphenate "Realworld" word if it is not already
% recTypes = regexprep(lower(recTypes),'realworld','real-world');

% Force capitalisation
% recTypes = regexprep(lower(recTypes),'(\<[a-z])','${upper($1)}');


colours = {[ ...            R G B  values
    1.0 0.0 0.0       ; ...       Active Wall Off
    0.2 0.2 1.0       ;]}; %      Active Wall On

lineStyles = { ...         Lines
    '--'       ; ...       Active Wall Off
    '- '       ;}; %       Active Wall On

markers = {[ ...          Marker Shapes
    '.'       ; ...       Active Wall Off
    '.'       ;]}; %      Active Wall On

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
% RT = numel(recTypes);
SN = numel(sigNames);
% Res_Matrix  = cell(RT,SN);
% Res_CI      = cell(RT,SN);
    
    
    %% Read results
    results_func = str2func(['Results.import_' 'RoomImpulseResponse_vs_Time']);
    
    % Read results
    [Xvec,Vals] ...
        = results_func( SYS );
    if Xvec == -1
        Results.Axes_Builders.Helpers.setErrorAxis(ax,Vals{3})
        return
    end
    
        for sn = 1:SN
            if any(sn==1:size(Vals,2))
                % Generate Plottable data matrices and vectors
                Res_Matrix{sn} = Vals(:,sn);
                Hrz_Vec = Xvec.';
                legendStrings = {legendStrings{:}, [sigNames{sn}]};
            else
                Res_Matrix{sn} = NaN(1,numel(Hrz_Vec));
                legendStrings = {legendStrings{:}, ['Not Found']};
            end
        end
    Res_Matrix = Res_Matrix.';
    
%%
cols = colours{:};
lines= lineStyles(:);
mrks = markers{:};

fs = SYS.signal_info.Fs;
[pkV,pkI] = max(Res_Matrix{1});
[pkVm,pkIm] = min(Res_Matrix{1});
domain = round(Hrz_Vec([0 0.01*fs] + pkI-floor(0.001*fs)),SYS.publication_info.sigRounding,'significant');
domain_lbl = 'Time ($\mathrm{s}$)';

axCurr = ax(1);
range = ceil(abs([pkVm pkV])*100)/100.*[-1 1];
if isfield(SYS.publication_info,'YLim_override')
    range = SYS.publication_info.YLim_override;
end
range_lbl = 'Normalised Amplitude';



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
    
%     ha=area([Hrz_Vec; Hrz_Vec].', [Res_Matrix{pl}.'+Res_CI{pl}(:,2), -sum(Res_CI{pl},2)], ...
%         PlDetailsShade{:} );
%     ha(1).FaceColor=cols(pl,:);drawnow;
%     ha(2).FaceColor=cols(pl,:);drawnow;
%     v=version;i=find(v=='.');V=str2num(v(1:(i(2)-1)));
%     if V >= 8.6 % If MATLAB version is greater than 8.6 (R2015b)
%         ha(1).FaceAlpha=0.0;
%         ha(2).FaceAlpha=trendAlpha;
%     else
%        ha(1).Face.ColorType = 'truecoloralpha';
%        ha(1).Face.ColorData(4) = 255 * 0.0;
%        ha(2).Face.ColorType = 'truecoloralpha';
%        ha(2).Face.ColorData(4) = 255 * trendAlpha;
%     end        
    
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

%% Some helpful summary values for publication
% disp('Simulated')
% disp(['Avg Quiet SPL: ' ...
%     num2str(mean(Res_Matrix{2}(Hrz_Vec>=SYS.signal_info.f_low/1e3 & Hrz_Vec<=1.601359980101278)))])
% disp(['Avg AC: ' ...
%     num2str(mean(Res_Matrix{1}(Hrz_Vec>=SYS.signal_info.f_low/1e3 & Hrz_Vec<=1.601359980101278)-Res_Matrix{2}(Hrz_Vec>=SYS.signal_info.f_low/1e3 & Hrz_Vec<=1.601359980101278))*2)])
% disp('Real-world')
% disp(['Avg Quiet SPL: ' ...
%     num2str(mean(Res_Matrix{4}(Hrz_Vec>=SYS.signal_info.f_low/1e3 & Hrz_Vec<=1.601359980101278)))])
% disp(['Avg AC: ' ...
%     num2str(mean(Res_Matrix{3}(Hrz_Vec>=SYS.signal_info.f_low/1e3 & Hrz_Vec<=1.601359980101278)-Res_Matrix{4}(Hrz_Vec>=SYS.signal_info.f_low/1e3 & Hrz_Vec<=1.601359980101278))*2)])

0; %breakpoint location

end



