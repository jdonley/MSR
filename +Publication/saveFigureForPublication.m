function saveFigureForPublication( SYS, figH )
%SAVEFIGUREFORPUBLICATION Summary of this function goes here
%
% Syntax:	SAVEFIGUREFORPUBLICATION( SYS, figH )
%
% Inputs:
% 	SYS - Soundfield reproduction system specification object
%   figH - Figure handle
%
% See also:

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016
% Date: 10 June 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DP = SYS.publication_info.DocumentPath;
FN = SYS.publication_info.FigureName;
FMT = SYS.publication_info.print_fmt;
RES = SYS.publication_info.print_res;

if ~exist(DP,'dir'); mkdir(DP); end
export_fig([DP '\' FN '.' FMT], ['-' FMT], ['-r' num2str(RES)], figH);

end

