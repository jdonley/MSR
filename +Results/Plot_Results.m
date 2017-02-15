function Plot_Results( SYS )
if nargin < 1
    % Load System
    SYS = Current_Systems.loadCurrentSRsystem;
end
close all;

%% Build Figure
figFunc = str2func( ['Results.Figure_Builders.' SYS.publication_info.FigureName] );
figH = figFunc( SYS );

%% Save Figure
Publication.saveFigureForPublication( SYS, figH );