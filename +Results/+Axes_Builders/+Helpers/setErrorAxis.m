function setErrorAxis(ax,errMSG,filepath)
%SETERRORAXIS Sets parameters of an axis from a SR_System object and other inputs
%
% Syntax:	SETERRORAXIS( ax, errMSG, filepath )
%
% Inputs:
%   ax - Axis handle to reuse
%   errMSG - Error message to display in the axis
%   filepath - Filepath to show in command window warning
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

axes(ax(1));
text(1,1, errMSG, 'HorizontalAlignment','center');

for i =1:numel(ax)
    axC = ax(i);
    axes(axC);
    axC.XTickMode = 'manual';
    axC.YTickMode = 'manual';
    axC.XTick = [];
    axC.YTick = [];
    axC.YColor = [0,0,0];
end
wrnCol = [255,100,0]/255;
cprintf(wrnCol, ['Warning: ' errMSG '\n''' strrep(filepath,'\','\\') '''\n']);

end