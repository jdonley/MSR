
clear;
close all;


%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

%% Build Figure
figFunc = str2func( ['Results.Figure_Builders.' SYS.publication_info.FigureName] );
figH = figFunc( SYS );

%% Save Figure
Publication.saveFigureForPublication( SYS, figH );