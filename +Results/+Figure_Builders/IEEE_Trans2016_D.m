function [ fig ] = IEEE_Trans2016_D( SYS )
%TEMPLATE Summary of this function goes here
%
% Syntax:	[OUTPUTARGS] = TEMPLATE(INPUTARGS) Explain usage here
%
% Inputs:
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
%
% Outputs:
% 	output1 - Description
% 	output2 - Description
%
% Example:
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
%
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016-2017
% Date: 28 November 2016
% Revision: 0.2 (23 March 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Name',SYS.publication_info.FigureName,'NumberTitle','off');
fig.Color = [1 1 1];

nCol = SYS.publication_info.subPlotDims(1);
nRow = SYS.publication_info.subPlotDims(2);
if numel(SYS.publication_info.subPlotDims)==3 ...
        && ~SYS.publication_info.MergeLines
    nPag = SYS.publication_info.subPlotDims(3);
else
    nPag = 1;
end

axHndls = tightPlots( ...
    nRow, nCol, ...
    SYS.publication_info.figure_width, ...
    SYS.publication_info.axis_aspect_ratio, ...
    SYS.publication_info.axes_gap, ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');

newAxs  = [];
legStrs = [];

%%
SYStmp = SYS;
for row = 1:nRow
    for col = 1:nCol
        for pag = 1:nPag
            
            I = sub2ind([nCol nRow],col,row);
            %             Ir  = sub2ind([nRow nCol],row,col);
            if SYS.publication_info.MergeLines
                nPag_ = SYS.publication_info.subPlotDims(3);
            else
                nPag_ = nPag;
            end
            Ir = arrayfun(@(x) sub2ind([nRow nCol nPag_],row,col,x), 1:nPag_);
            
            SYStmp.Main_Setup   = SYS.Main_Setup(Ir);
            SYStmp.Masker_Setup = SYS.Masker_Setup(Ir);
            
            SYStmp.signal_info.method = ...
                SYS.signal_info.methods_list{ ...
                SYS.signal_info.methods_list_masker(Ir)};
            
            % Set a tag to distinguish between axes
            axes(axHndls( I )); ax = gca; ax.Tag = char(64+I);
            
            [nA, lS] = Results.Axes_Builders.STOI_PESQ( ...
                SYStmp, ...
                axHndls( I ) );
            
        end
        newAxs  = [newAxs; nA];
        legStrs = [legStrs; lS];
        
    end
    
end
SYS_ = SYS;
SYS_.Main_Setup(row*col+1:end)=[];
SYS_.Masker_Setup(row*col+1:end)=[];

drawnow;pause(0.05);
tightfig(fig);

Results.Figure_Builders.Helpers.setLegend( SYS_, newAxs, legStrs, 'errorbar' );

Results.Figure_Builders.Helpers.setLabels( SYS_, newAxs, nRow, nCol );

end



