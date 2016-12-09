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

ax.XLim = domain + [-1 1]*diff(domain)*axBuf(1);
ax.YLim = range + [-1 1]*diff(range)*axBuf(2);

ax.TickDir = SYS.publication_info.axes_tickdir; % Set tick direction

ax.XTickMode = 'manual';
ax.YTickMode = 'manual';

ax.XTick = round( domain(1):diff(domain)/(nXT-1):domain(end) , ...
    sigRounding, 'significant');
ax.YTick = round( range(1):diff(range)/(nYT-1):range(end), ...
    sigRounding, 'significant');

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

if strcmp(SYS.publication_info.axes_grid, 'minor') %if minor then turn on major and minor
    grid(ax, 'on');
end
grid(ax, SYS.publication_info.axes_grid);

end