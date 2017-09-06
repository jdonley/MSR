function [ fig ] = IEEE_Letters2017_A( SYS )
%TEMPLATE Summary of this function goes here
%
% Syntax:	[ fig ] = IEEE_Letters2017_A( SYS )
%
% Inputs:
% 	SYS - Soundfield Reproduction system object
%
% Outputs:
% 	fig - The publication ready figure handle
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
% Copyright: Jacob Donley 2017
% Date: 4 September 2017
% Version: 0.1 (4 September 2017)
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
        
        if numel(SYS.Room_Setup) > 1
            % If there is more than one room then we choose the room set up for
            % reproduction (transmission)
            SYS.Room_Setup = SYS.Room_Setup( ...
                strcmpi({SYS.Room_Setup.SystemType},'transmit'));
        end
        
        axNo = sub2ind([Nrow,Ncol],c,r);
        [nA, lS] = Results.Axes_Builders.SUPPRESSION_Freq( ...
            SYS, ...
            axHndls( axNo ) );
        newAxs  = [newAxs; nA];
        legStrs = [legStrs; lS];
        
    end
    
end


%%
drawnow;pause(0.05);
tightfig(fig);

Results.Figure_Builders.Helpers.setLegend( SYS, newAxs, legStrs );

% Results.Figure_Builders.Helpers.setLabels( SYS, newAxs, Nrow, Ncol );


end



