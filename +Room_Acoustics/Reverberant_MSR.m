clear;
%clear classes;
%close all;
fclose all;clc;
delete(gcp('nocreate'));

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

%%
for typ = 1:length(SYS.signal_info.methods_list)
    
    SYS.signal_info.method = SYS.signal_info.methods_list{typ};
    
    Room_Acoustics.Apply_RIRs.Reverberant_MSR_batchfunc( SYS );
    
end