function setLegend( SYS, axs, strs, linetypes )
%SETLEGEND Sets up the legend in the current figure
%
% Syntax:	SETLEGEND( SYS, axs, strs, linetypes )
%
% Inputs:
% 	SYS - Soundfield reproduction system specification object
% 	axs - Axis handles to obtain legend information from
% 	strs - Legend entry strings
% 	linetypes - Line types to show in the legend
%
%
% See also: 

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016
% Date: 14 June 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legEntrySpacing = 0.0; % percentage of entry width

if nargin < 4
    linetypes = 'Line';
end
if diff(size(axs))<1
    axs = axs.';
end
axes( axs(1) ); % Pick the first axes
ax = gca;

nLegEnts = size(strs,2);
if ~nLegEnts, return; end

% Search all axes
axEnts=[];
for a = 1:numel(axs)
    axEnts = [axEnts; findobj(axs(a).Children,'Type',linetypes)];
end
% Just keep the first few entries
axEnts = axEnts(1:nLegEnts);

FF = SYS.publication_info.LaTeX_FontFamily; % Set fonts for latex interpreter
FS = SYS.publication_info.FontSize;
latexFontSettings = ['{\fontfamily{' FF '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];
for l = 1:nLegEnts
    if strcmpi(SYS.publication_info.Interpreter, 'latex')
        strs{1,l} = [latexFontSettings ...
            strs{1,l} '}'];
    end
end

% Create legend
[leg,legi] = legend(ax, axEnts ,strs(1,:), ...
    'Interpreter', SYS.publication_info.Interpreter, ...
    'FontName', SYS.publication_info.FontName, ...
    'FontSize', SYS.publication_info.FontSize);
% [leg] = legend(ax, axEnts ,strs(1,:), ...
%     'Interpreter', SYS.publication_info.Interpreter, ...
%     'FontName', SYS.publication_info.FontName, ...
%     'FontSize', SYS.publication_info.FontSize);
leg.Box = 'off';
legli = findobj(legi,'Type','line'); % legend lines

% Change legend grid
tmpUnits = leg.Units; leg.Units = 'points';
legWid = leg.Position(3) * (1-diff(legli(2).XData)) ... %normalised units
    + max([legli(1:2:end).MarkerSize]); % plus maximum marker width in points
legWid = legWid * (1+legEntrySpacing); % Add spacing between legend entries as percentage of legWid
gridLegend(axEnts, nLegEnts, strs(1,:), ...
    'existinglegend',leg, ...
    'existinglegendicons',legi, ...
    'location','stationary', ...
    'legendwidth',legWid);

% Resize legend markers
for ll = 1:2:length(legli)
legli(ll).XData(1)  = legli(ll).XData(1) ...
    + diff([legli(ll).XData(1), legli(ll+1).XData(2)]) .* 1;
end

% Find max width and max height of axes
allAxs = findobj(ax.Parent.Children,'type','axes');
axPos=reshape([allAxs.Position],4,[])';
fullWidHigh = max(axPos(:,1:2)) + max(axPos(:,3:4));

% Relocate legend
% leg.Position(1:2) = fullWidHigh - [leg.Position(3) 0]; % upper right corner
leg.Position(1:2) = fullWidHigh .* [1/2 1] - leg.Position(3:4).*[1/2 1/nLegEnts]; % centered below title
leg.Units = tmpUnits;

% Resize entire figure
titTI_orig = ax.TightInset;
nXtraLines =  nLegEnts-1;
if isempty(ax.Title.String),ax.Title.String={ax.Title.String};end
ax.Title.String = [ax.Title.String{1}; ...
    mat2cell( repmat( ' ', nXtraLines, 1 ), ones(nXtraLines,1) );] ; % minus 1 is for the title
ax.Parent.Position(4) = ax.Parent.Position(4) + diff( [titTI_orig(4), ax.TightInset(4)] );

end