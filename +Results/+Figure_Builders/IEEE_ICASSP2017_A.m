function [ fig ] = IEEE_ICASSP2017_A( SYS )
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

Nrow = numel(SYS.signal_info.recording_type);
Ncol = size(SYS.signal_info.methods_list_masker,2);

axHndls = tightPlots( ...
    Nrow, Ncol, ...
    SYS.publication_info.figure_width, ...
    SYS.publication_info.axis_aspect_ratio, ...
    SYS.publication_info.axes_gap, ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');

newAxs  = [];
legStrs = [];

%%

for r = 1:Nrow
    
    for c = 1:Ncol
        
        SYS.signal_info.method = SYS.signal_info.methods_list{end};
        
        axNo = sub2ind([Nrow,Ncol],c,r);
        [nA, lS] = Results.Axes_Builders.SUPPRESSION( ...
            SYS, ...
            axHndls( axNo ) );
        newAxs  = [newAxs; nA];
        legStrs = [legStrs; lS];
        
    end
    
end


%%
drawnow;pause(0.05);
tightfig(fig);

Results.Figure_Builders.Helpers.setLegend( SYS, newAxs, legStrs, 'errorbar' );

% Results.Figure_Builders.Helpers.setLabels( SYS, newAxs, Nrow, Ncol );


end



