function [ fig ] = IEEE_Trans2016_C( SYS )
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
% Copyright: Jacob Donley 2016
% Date: 09 June 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Name',SYS.publication_info.FigureName,'NumberTitle','off');
fig.Color = [1 1 1];

nCol = SYS.publication_info.subPlotDims(1);
nRow = SYS.publication_info.subPlotDims(2);

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
    SYStmp.signal_info.recording_type = ...
        SYS.signal_info.recording_type{row};
    
    for col = 1:nCol
        I = sub2ind([nCol nRow],col,row);
        
        SYStmp.Main_Setup   = SYS.Main_Setup( I );
        SYStmp.Masker_Setup = SYS.Masker_Setup( I );
        
        SYStmp.signal_info.method = ...
            SYS.signal_info.methods_list{ ...
            SYS.signal_info.methods_list_masker( I )};
        
        % Set a tag to distinguish between axes
        axes(axHndls( I )); ax = gca; ax.Tag = char(64+I);
        
        [nA, lS] = Results.Axes_Builders.STOI_PESQ( ...
            SYStmp, ...
            axHndls( I ) );
        
        newAxs  = [newAxs; nA];
        legStrs = [legStrs; lS];
    end
    
end

drawnow;pause(0.05);
tightfig(fig);

Results.Figure_Builders.Helpers.setLegend( SYS, newAxs, legStrs, 'errorbar' );

Results.Figure_Builders.Helpers.setLabels( SYS, newAxs, nRow, nCol );

end



