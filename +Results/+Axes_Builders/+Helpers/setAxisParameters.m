function setAxisParameters( SYS, ax, range, domain, range_lbl, domain_lbl)
%SETAXISPARAMETERS Sets parameters of an axis from a SR_System object and other inputs
%
% Syntax:	SETAXISPARAMETERS( SYS, ax, range, domain, range_lbl, domain_lbl )
%
% Inputs:
% 	SYS - Soundfield reproduction system specification object
%   ax - Axis handle to reuse
%   range - Axis range (two element vector)
%   domain - Axis domain (two element vector)
%   range_lbl - Axis y-label
%   domain_lbl - Axis x-label
%
%
% See also:

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016
% Date: 23 November 2016
% Revision: 0.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axBuf = SYS.publication_info.axes_limitBufs; % axis limits buffer in percentage [width, height]
nXT = SYS.publication_info.axes_NumTicks(1); % number of X ticks
nYT = SYS.publication_info.axes_NumTicks(2); % number of X ticks
sigRounding = SYS.publication_info.sigRounding; % number of significant figures rounding
FF = SYS.publication_info.LaTeX_FontFamily; % Set fonts for latex interpreter
FFnums = SYS.publication_info.LaTeX_NumbersFontFamily; % Set number fonts for latex interpreter
FS = SYS.publication_info.FontSize;

latexFontSettings = ['{\fontfamily{' FF '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];
latexNumFontSettings = ['{\fontfamily{' FFnums '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];

ax.FontName = SYS.publication_info.NumbersFontName;
ax.FontSize = SYS.publication_info.FontSize;

ax.Title.Interpreter = SYS.publication_info.Interpreter;
ax.Title.FontName = SYS.publication_info.FontName;
ax.Title.FontSize = SYS.publication_info.FontSize;
if strcmpi(ax.Title.Interpreter, 'latex')
    titleTxt = [latexFontSettings ...
        SYS.publication_info.FigureTitle '}'];
else
    titleTxt =  SYS.publication_info.FigureTitle;
end
ax.Title.String = {titleTxt; ...
    ' '}; % Extra line is for axes letters

ax.XLabel.Interpreter = SYS.publication_info.Interpreter;
ax.YLabel.Interpreter = SYS.publication_info.Interpreter;
ax.XLabel.FontName = SYS.publication_info.FontName;
ax.YLabel.FontName = SYS.publication_info.FontName;
ax.XLabel.FontSize = SYS.publication_info.FontSize;
ax.YLabel.FontSize = SYS.publication_info.FontSize;

ax.XLabel.String = domain_lbl;
ax.YLabel.String = range_lbl;

if strcmpi(ax.XLabel.Interpreter, 'latex')
    ax.XLabel.String = [latexFontSettings ax.XLabel.String '}'];
end
if strcmpi(ax.YLabel.Interpreter, 'latex')
    ax.YLabel.String = [latexFontSettings ax.YLabel.String '}'];
end

if isfield(SYS.publication_info,'axes_Scales')
    ax.XScale = SYS.publication_info.axes_Scales{1};
    ax.YScale = SYS.publication_info.axes_Scales{2};    
end

ax.TickDir = SYS.publication_info.axes_tickdir; % Set tick direction

ax.XTickMode = 'manual';
ax.YTickMode = 'manual';

if ~isempty(strfind(ax.XScale, 'lin'))
    XTick = round( domain(1):diff(domain)/(nXT-1):domain(end) , ...
        sigRounding, 'significant');
    if any(diff(XTick)==0)
        XTick = domain(1):diff(domain)/(nXT-1):domain(end);
    end
    ax.XTick = XTick;
elseif ~isempty(strfind(ax.XScale, 'log'))
    XTick = 10.^round( log10(domain(1)):diff(log10(domain))/(nXT-1):log10(domain(end)) , ...
        sigRounding, 'significant');
    if any(diff(XTick)==0)
        XTick = 10.^ ( log10(domain(1)):diff(log10(domain))/(nXT-1):log10(domain(end)) );
    end
    ax.XTick = XTick;
end
if ~isempty(strfind(ax.YScale, 'lin'))
    ax.YTick = round( range(1):diff(range)/(nYT-1):range(end), ...
        sigRounding, 'significant');
elseif ~isempty(strfind(ax.YScale, 'log'))
    ax.YTick = 10.^round( log10(range(1)):diff(log10(range))/(nYT-1):log10(range(end)), ...
        sigRounding, 'significant');
end

if isfield(SYS.publication_info,'XTicks_override')
    ax.XTick = SYS.publication_info.XTicks_override;
    domain = SYS.publication_info.XTicks_override([1 end]);
end
if isfield(SYS.publication_info,'YTicks_override')
    ax.YTick = SYS.publication_info.YTicks_override;
    range = SYS.publication_info.YTicks_override([1 end]);
end

if ~isempty(strfind(ax.XScale, 'lin'))
    diffSepDomain = diff(domain)*axBuf(1);
elseif ~isempty(strfind(ax.XScale, 'log'))
    diffSepDomain = diff(domain)*axBuf(1);
end
if ~isempty(strfind(ax.YScale, 'lin'))
    diffSepRange = diff(range)*axBuf(2);
elseif ~isempty(strfind(ax.YScale, 'log'))
    diffSepRange = diff(range)*axBuf(2);
end

ax.XLim = domain + [-1 1]*diffSepDomain;
ax.YLim = range + [-1 1]*diffSepRange;

ax.TickLabelInterpreter = SYS.publication_info.Interpreter;
if strcmpi(ax.TickLabelInterpreter, 'latex')
    for t = 1:numel(ax.XTickLabel)
        ax.XTickLabel(t) = {[latexNumFontSettings ax.XTickLabel{t} '}']};
    end
    for t = 1:numel(ax.YTickLabel)
        ax.YTickLabel(t) = {[latexNumFontSettings ax.YTickLabel{t} '}']};
    end
end

ax.YColor = [0,0,0]; %Hard coded black axes colors

grid(ax, 'off');
if strcmp(SYS.publication_info.axes_grid, 'minor') %if minor then turn on major and minor
    grid(ax, 'on');
end
grid(ax, SYS.publication_info.axes_grid);
ax.YMinorGrid = SYS.publication_info.axes_gridMinorY;
ax.XMinorGrid = SYS.publication_info.axes_gridMinorX;

end