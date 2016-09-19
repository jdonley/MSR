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

nRT = numel(SYS.signal_info.recording_type);
nML = size(SYS.signal_info.methods_list_masker,2);

axHndls = tightPlots( ...
    nRT, nML, ...
    SYS.publication_info.figure_width, ...
    SYS.publication_info.axis_aspect_ratio, ...
    SYS.publication_info.axes_gap, ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');

newAxs  = [];
legStrs = [];

%%
RTs_tmp = SYS.signal_info.recording_type;
for rt = 1:nRT
    SYS.signal_info.recording_type = ...
        RTs_tmp{rt};
    
    for ml = 1:nML
        SYS.signal_info.method = ...
            SYS.signal_info.methods_list{ ...
            SYS.signal_info.methods_list_masker(:,ml)};
        
        axNo = sub2ind([nRT,nML],ml,rt);
        [nA, lS] = Results.Axes_Builders.STOI_PESQ( ...
            SYS, ...
            axHndls( axNo ) );
        newAxs  = [newAxs; nA];
        legStrs = [legStrs; lS];
        
    end
    
end
SYS.signal_info.recording_type = RTs_tmp;

drawnow;pause(0.05);
tightfig(fig);

Results.Figure_Builders.Helpers.setLegend( SYS, newAxs, legStrs, 'errorbar' );

Results.Figure_Builders.Helpers.setLabels( SYS, newAxs, nRT, nML );

end



