clear;
%clear classes;
%close all;
fclose all;
delete(gcp('nocreate'));

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

%%
for typ = 1:length(SYS.signal_info.methods_list)
    
    SYS.signal_info.method = SYS.signal_info.methods_list{typ};
    
    subSYS = SYS;
    if any(typ == subSYS.signal_info.methods_list_clean)
        subSYS.Main_Setup = subSYS.Main_Setup(typ);
    elseif any(typ == subSYS.signal_info.methods_list_masker)
        subSYS.Main_Setup = subSYS.Masker_Setup(typ);
    end
    
    Room_Acoustics.Apply_RIRs.Reverberant_MSR_batchfunc( subSYS );
    
end
