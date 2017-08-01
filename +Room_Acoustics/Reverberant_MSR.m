function Reverberant_MSR(SYS)
% clear;
%clear classes;
%close all;
fclose all;
delete(gcp('nocreate'));

%% Load System
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end

%%

N = length(SYS.signal_info.methods_list_clean(SYS.signal_info.methods_list_clean>=1)) ...
    + length(SYS.signal_info.methods_list_masker(SYS.signal_info.methods_list_masker>=1));
paired = isfield(SYS.signal_info,'methods_list_paired') && SYS.signal_info.methods_list_paired;

for typ = [2 6]%1:N
    
    SYS.signal_info.method = SYS.signal_info.methods_list{typ};
    
    subSYS = SYS;
    if any(typ == subSYS.signal_info.methods_list_clean)
        subSYS.Main_Setup = subSYS.Main_Setup(typ == subSYS.signal_info.methods_list_clean);
        if paired
            subSYS.Masker_Setup = subSYS.Masker_Setup(typ == subSYS.signal_info.methods_list_clean);
        end
    elseif any(typ == subSYS.signal_info.methods_list_masker)
        subSYS.Masker_Setup = subSYS.Masker_Setup(typ == subSYS.signal_info.methods_list_masker);
        if paired
            subSYS.Main_Setup = subSYS.Main_Setup(typ == subSYS.signal_info.methods_list_masker);
        end
    end
    
    if isfield(SYS.signal_info, 'UseMeasuredATFs') && SYS.signal_info.UseMeasuredATFs, Room_Acoustics.useMeasuredATF(subSYS); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Room_Acoustics.Apply_RIRs.Reverberant_MSR_batchfunc( subSYS );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isfield(SYS.signal_info, 'UseMeasuredATFs') && SYS.signal_info.UseMeasuredATFs, Room_Acoustics.useSimulatedATF(subSYS); end % Revert back after use
    
end


