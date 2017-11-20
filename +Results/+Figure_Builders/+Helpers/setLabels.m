function setLabels( SYS, axs, nY, nX )
%SETLABELS Sets up the axes labels in the current figure
%
% Syntax:	SETLABELS( SYS, axs, nY, nX )
%
% Inputs:
% 	SYS - Soundfield reproduction system specification object
% 	axs - Axis handles to obtain legend information from
% 	nY - Number of axes in Y direction
% 	nX - Number of axes in X direction
%
%
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016
% Date: 14 June 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FF = SYS.publication_info.LaTeX_FontFamily; % Set fonts for latex interpreter
FS = SYS.publication_info.FontSize;
latexFontSettings = ['{\fontfamily{' FF '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];

LblAxInds = []; %Indices of axes with labels [title; y1Lbl; xLbl; y2Lbl;]
for y = 1:nY
    for x = 1:nX
        for lr = 1:size(axs,2) %left right yaxis
            axInd = [sub2ind([nX,nY],x,y), lr];
            axes( axs( axInd(1), axInd(2) ) );
            ax = gca;
            
            
            if y ~= nY % if not the bottom row
                ax.XTickLabel = '';
            end
            if x ~= 1 && lr == 1 % if not the left column and left y axis
                ax.YTickLabel = '';
            end
            if x ~= nX && lr == 2 % if not the right column and right y axis
                ax.YTickLabel = '';
            end
            
            if lr == size(axs,2) %Axes numbers
                if isfield(SYS.publication_info,'subPlotTitles')
                    addTitleTxt = [': ' SYS.publication_info.subPlotTitles{axInd(1)}];
                else 
                    addTitleTxt = [];
                end
                ax.Title.String = ['(' char( 64 + axInd(1) ) ')' addTitleTxt];

                if strcmpi(SYS.publication_info.Interpreter, 'latex')
                    ax.Title.Interpreter = SYS.publication_info.Interpreter;
                    ax.Title.String = [latexFontSettings ...
                        ax.Title.String '}'];
                end
                tmpUnits = ax.Title.Units; ax.Title.Units = 'points'; %assume points
                ax.Title.HorizontalAlignment = 'left';
                ax.Title.Position(1) = 0; 
                if (contains(lower(SYS.publication_info.axes_tickdir),'both') ...
                        || contains(lower(SYS.publication_info.axes_tickdir),'out'))
                    ax.Title.Position(2) = ax.Position(4);
                end
                ax.Title.Units = tmpUnits; % 0.05 => 5 percent of axes width
                %ax.Title.Position(1) = ax.Position(3)*0.05; ax.Title.Units = tmpUnits; % 0.05 => 5 percent of axes width
            end
            if ~(y == 1 && x == 1) && lr == 1 % if not top left axis and left y axis
                ax.Title.String = '';
            elseif y == 1 && x == 1 && lr == 1
                LblAxInds(1,:) = axInd;
            end
            
            if ~(y == nY && x == 1) && lr == 1 % if not bottom left axis and left y axis
                ax.XLabel.String = '';
                ax.YLabel.String = '';
            elseif y == nY && x == 1 && lr == 1
                LblAxInds(2,:) = axInd;
                LblAxInds(3,:) = axInd;
            end
            if ~(y == nY && x == nX) && lr == 2 % if not bottom right axis and right y axis
                ax.YLabel.String = '';
            elseif y == nY && x == nX && lr == 2
                LblAxInds(4,:) = axInd;
            end
        end
    end
end

allAxs = findobj(ax.Parent.Children,'type','axes');
axPos=reshape([allAxs.Position],4,[])';
fullWidHigh = diff([min(axPos(:,1:2)); ...
    (max(axPos(:,1:2)) + max(axPos(:,3:4))) ]);

% Centre each axis label
for lbl = 1:size(LblAxInds,1)
    axes( axs(LblAxInds(lbl,1),LblAxInds(lbl,2)) );
    ax = gca;
    tmpObj = [];
    if mod(lbl,2)
        if lbl == 1
            tmpObj = ax.Title;
        else
            tmpObj = ax.XLabel;
        end
    else
        tmpObj = ax.YLabel;
    end
    
    tmpUnit = tmpObj.Units; tmpObj.Units = 'points'; % assume units are points
    tmpObj.Position(2-mod(lbl,2)) = fullWidHigh(2-mod(lbl,2))/2;
    tmpObj.Units = tmpUnit;
    ax.LooseInset = ax.TightInset;
    
    if size(axs,2) > 1 && LblAxInds(lbl,2) == 1 %if a left y axis then reset right y axis on top
        axes( axs(LblAxInds(lbl,1), 2 ) );
    end
end

end