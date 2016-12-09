clear;
%clear classes;
%close all;
fclose all;
delete(gcp('nocreate'));

%% Load System
SYS = Current_Systems.loadCurrentSRsystem;

%%
paired = isfield(SYS.signal_info,'methods_list_paired') && SYS.signal_info.methods_list_paired;

ml_tmp = SYS.signal_info.methods_list;
rt_tmp = SYS.signal_info.recording_type;
Nrt = numel(SYS.signal_info.recording_type);
for rt = 1:Nrt
    SYS.signal_info.recording_type = rt_tmp{rt};

    for c = SYS.signal_info.methods_list_clean
        masker_list = SYS.signal_info.methods_list_masker;
        if paired
            masker_list = masker_list(c);
        end
        for m = masker_list
            m(m<1)=[];
            
            SYS.signal_info.methods_list = {...
                ml_tmp{c}, ...
                ml_tmp{m} };
            
            subSYS = SYS;
            subSYS.Main_Setup(~(c==SYS.signal_info.methods_list_clean))=[];
            subSYS.Masker_Setup(~(m==SYS.signal_info.methods_list_masker))=[];
            if paired
                subSYS.signal_info.methods_list_clean = 1:numel(c);
                subSYS.signal_info.methods_list_masker = numel(c)+1:numel(c)+numel(m);
            end
            
            Room_Acoustics.Apply_RIRs.Reverberant_MSR_Analysis_batchfunc( subSYS );
            
        end
    end
end
SYS.signal_info.recording_type = rt_tmp;
SYS.signal_info.methods_list = ml_tmp;
    
