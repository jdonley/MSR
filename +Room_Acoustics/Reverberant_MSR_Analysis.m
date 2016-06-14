clear;
%clear classes;
%close all;
fclose all;clc;
delete(gcp('nocreate'));

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

%%
ml_tmp = SYS.signal_info.methods_list;
rt_tmp = SYS.signal_info.recording_type;
for rt = 1:numel(SYS.signal_info.recording_type)
    SYS.signal_info.recording_type = rt_tmp{rt};

    for c = SYS.signal_info.methods_list_clean
        for m = SYS.signal_info.methods_list_masker
            m(m<1)=[];
            
            SYS.signal_info.methods_list = {...
                ml_tmp{c}, ...
                ml_tmp{m} };
            Room_Acoustics.Apply_RIRs.Reverberant_MSR_Analysis_batchfunc( SYS );
        
        end
    end
end
SYS.signal_info.recording_type = rt_tmp;
SYS.signal_info.methods_list = ml_tmp;
    
