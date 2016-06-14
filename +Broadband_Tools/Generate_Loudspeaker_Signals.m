clear;
%clear classes;
close all;
fclose all;clc;

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

%%

for typ = 1:length(SYS.signal_info.methods_list)
    
    SYS.signal_info.method = SYS.signal_info.methods_list{typ};
    
    Broadband_Tools.Generate_Loudspeaker_Signals_batchfunc( SYS );
    
end