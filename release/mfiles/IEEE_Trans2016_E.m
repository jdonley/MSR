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
% Copyright: Jacob Donley 2016
% Date: 28 November 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fig = figure('Name',SYS.publication_info.FigureName,'NumberTitle','off');
fig.Color = [1 1 1];

nRow = numel(SYS.Main_Setup);
nCol = size(SYS.signal_info.methods_list_masker,2);

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
                
        I = sub2ind([nCol nRow],col,row);

        SYStmp.Main_Setup   = SYS.Main_Setup( I );
        
        SYStmp.signal_info.method = SYS.signal_info.methods_list{end};
%         SYS.signal_info.recording_type = SYS.signal_info.recording_type{r};
                
        [nA, lS] = Results.Axes_Builders.SPL_Freq( ...
            SYStmp, ...
            axHndls( I ) );
        newAxs  = [newAxs; nA];
        legStrs = [legStrs; lS];
        
    end
    
end


%%
drawnow;pause(0.05);
tightfig(fig);

Results.Figure_Builders.Helpers.setLegend( SYS, newAxs, legStrs );

Results.Figure_Builders.Helpers.setLabels( SYS, newAxs, nRow, nCol );


end
