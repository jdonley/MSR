function Reverberant_MSR_Analysis(SYS)
% clear;
%clear classes;
%close all;
fclose all;
delete(gcp('nocreate'));

%% Load System
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end

%%
paired = isfield(SYS.signal_info,'methods_list_paired') && SYS.signal_info.methods_list_paired;

ml_tmp = SYS.signal_info.methods_list;
rt_tmp = SYS.signal_info.recording_type;
Nrt = numel(SYS.signal_info.recording_type);
for rt = 1:Nrt
    SYS.signal_info.recording_type = rt_tmp{rt};

    for c_ = 1:numel(SYS.signal_info.methods_list_clean)
        c = SYS.signal_info.methods_list_clean(c_);
        if isfield(SYS.signal_info,'recording_type_list') ... <- This retains backwards compatibility
                && ~strcmpi(...
                    SYS.signal_info.recording_type_list(c), ...
                    SYS.signal_info.recording_type) ...
                && paired
                % If the paired setups are not of the correct recording 
                % type then skip the analysis
                continue; 
        end
        masker_list = SYS.signal_info.methods_list_masker;
        if paired
            masker_list = masker_list(c_);
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
            else
                subSYS.signal_info.methods_list_clean = 1;
            end
            if isfield(SYS.system_info,'CurrentSpeakerArrayType') ... <- This retains backwards compatibility
                    && ...
                    ~strcmpi( SYS.system_info.CurrentSpeakerArrayType, ...
                    subSYS.Main_Setup.Speaker_Array_Type) ...
                    && (SYS.signal_info.UseMeasuredATFs ...
                    || strcmpi(SYS.signal_info.recording_type,'real-world'))
                % If the objects array type is not the current physically
                % set up type in the real-world and we are using measured
                % ATFs (which use that current physical setup) or the 
                % recording type is 'real-world' (which then must require a 
                % different physical setup) then we skip the analysis to
                % save time.
                continue; 
            end
            
            Room_Acoustics.Apply_RIRs.Reverberant_MSR_Analysis_batchfunc( subSYS );
            
        end
    end
end
SYS.signal_info.recording_type = rt_tmp;
SYS.signal_info.methods_list = ml_tmp;
    
