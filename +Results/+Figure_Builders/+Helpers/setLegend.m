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

if isfield(SYS.publication_info, 'showlegend') && ~SYS.publication_info.showlegend
    return;
end

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
nAxs = numel(axs);
for a = 1:nAxs
    axEnts = [axEnts; findobj(axs(a).Children,'Type',linetypes)];
end
axEnts(strcmpi({axEnts.Tag},'noLegend')) = [];
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
legtx = findobj(legi,'Type','text'); % legend lines
if isfield(SYS.publication_info,'leg_MarkerSize')
    for i = 1:numel(legli)
        legli(i).MarkerSize = SYS.publication_info.leg_MarkerSize;
    end
end
% Center align legend text
maxTxtWid = max(reshape([legtx.Extent],4,[])',[],1);
maxTxtWid = maxTxtWid(3);
for i = 1:numel(legtx)
    legtx(i).HorizontalAlignment = 'center';    % Center the text
    legtx(i).Position(1) = legtx(i).Position(1) + maxTxtWid/2; % re-position text block to center
end


% Change legend grid
tmpUnits = leg.Units; leg.Units = 'points';
if strcmpi(linetypes,'line'), a = 1; b=0.65; else a = 2; b=0.0; end
legWid = leg.Position(3) * (1-diff(legli(a).XData)*(1-b)) ... %normalised units
    + max([legli(1:2:end).MarkerSize]); % plus maximum marker width in points
legWid = legWid * (1+legEntrySpacing); % Add spacing between legend entries as percentage of legWid
if legWid*nLegEnts > axEnts(1).Parent.Parent.OuterPosition(3)
    LegCols = floor((axEnts(1).Parent.Parent.OuterPosition(3) / legWid)/2)*2;
else
    LegCols = nLegEnts;
end
gridLegend(axEnts, LegCols, strs(1,:), ...
    'existinglegend',leg, ...
    'existinglegendicons',legi, ...
    'location','stationary', ...
    'legendwidth',legWid);

leg.Position(4) = leg.Position(4)/LegCols;

% Resize legend markers
if ~strcmpi(linetypes,'line')
    for ll = 1:2:length(legli)
%         legli(ll).XData(1) = legli(ll).XData(1) ...
%             + diff([legli(ll).XData(1), legli(ll+1).XData(2)]) .* 1;
        legli(ll).XData(1) = legtx((ll+1)/2).Extent(1);
        unitTmp = legtx((ll+1)/2).Units; legtx((ll+1)/2).Units = 'points';
        legtx((ll+1)/2).Position(1) = legtx((ll+1)/2).Position(1) ...
            + legli(ll).MarkerSize;
        legtx((ll+1)/2).Units = unitTmp;
    end
else
    for ll = 1:2:length(legli)
        legli(ll).XData = [-diff(legli(ll).XData).*b 0] + legtx((ll+1)/2).Extent(1);
        legli(ll+1).XData(1) = legli(ll).XData(1) ...
            + diff(legli(ll).XData) .* 0.8;
    end
end

% Find max width and max height of axes
allAxs = findobj(ax.Parent,'type','axes');
% axPos=reshape([allAxs.OuterPosition],4,[])';
% fullWidHigh = max(axPos(:,1:2),[],1) ...
%             + max(axPos(:,3:4),[],1);
% if strcmpi(linetypes,'line')
%     fullWidHigh(1)  = fullWidHigh(1) + max(axPos(:,1));
% end
axPos=reshape([allAxs.Position],4,[])';
fullWidHigh = max(axPos(:,1:2) + axPos(:,3:4),[],1);

% Relocate legend
% leg.Position(1:2) = fullWidHigh - [leg.Position(3) 0]; % upper right corner
leg.Position(1:2) = fullWidHigh .* [1/2 1.0] - leg.Position(3:4).*[1/2 0.5/LegCols]; % centered below title
tmpUnits2 = allAxs(end).Title.Units; allAxs(end).Title.Units = 'points';
leg.Position(2) = allAxs(end).Position(2) + allAxs(end).Title.Position(2) - leg.Position(4)/2 + FS/2;
allAxs(end).Title.Units = tmpUnits2;
leg.Units = tmpUnits;

% Resize entire figure
titTI_orig = ax.TightInset;
nXtraLines =  LegCols;
if isempty(ax.Title.String),ax.Title.String={ax.Title.String};end

tmpUnits={};
tmpUnits{1} = ax.Title.Units; ax.Title.Units='points';
tmpUnits{2} = leg.Units; leg.Units = 'points';

ax.Title.String(2:end)=[];
% if nLegEnts > 1
%  ax.Title.String = [ax.Title.String{1}; {' '};] ;
% end
ax.Title.Position(2) = axEnts(1).Parent.Position(4) + leg.Position(4) + FS * ((nAxs > 1));
leg.Position(2) = axEnts(1).Parent.Position(2) + axEnts(1).Parent.Position(4) + FS * ((nAxs > 1));
ax.Parent.Position(4) = ax.Position(2) + ax.Title.Position(2) + FS;

ax.Title.Units=tmpUnits{1}; leg.Units = tmpUnits{2};
end