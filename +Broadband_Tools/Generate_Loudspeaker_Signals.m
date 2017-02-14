function Generate_Loudspeaker_Signals(SYS)
% clear;
%clear classes;
% close all;
fclose all;

%% Load System
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end

%%
N = length(SYS.signal_info.methods_list_clean(SYS.signal_info.methods_list_clean>=1)) ...
    + length(SYS.signal_info.methods_list_masker(SYS.signal_info.methods_list_masker>=1));
paired = isfield(SYS.signal_info,'methods_list_paired') && SYS.signal_info.methods_list_paired;

for typ = 1:N
    
    SYS.signal_info.method = SYS.signal_info.methods_list{typ};
    
    subSYS = SYS;
    if any(typ == subSYS.signal_info.methods_list_clean)
        if ~all(size(subSYS.Main_Setup)==1)
            subSYS.Main_Setup = subSYS.Main_Setup(typ == subSYS.signal_info.methods_list_clean);
        else
            subSYS.Main_Setup = subSYS.Main_Setup;
        end
        if paired
            subSYS.Masker_Setup = subSYS.Masker_Setup(typ == subSYS.signal_info.methods_list_clean);
        end
    elseif any(typ == subSYS.signal_info.methods_list_masker)
        if ~all(size(subSYS.Masker_Setup)==1)
            subSYS.Masker_Setup = subSYS.Masker_Setup(typ == subSYS.signal_info.methods_list_masker);
        else
            subSYS.Masker_Setup = subSYS.Masker_Setup;
        end
        if paired
            subSYS.Main_Setup = subSYS.Main_Setup(typ == subSYS.signal_info.methods_list_masker);
        end
    end
    
    Broadband_Tools.Generate_Loudspeaker_Signals_batchfunc( subSYS );
    
end

