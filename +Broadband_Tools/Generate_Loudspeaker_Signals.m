function Generate_Loudspeaker_Signals(SYS)
% Summary of this function goes here
% 
% Syntax:	[OUTPUTARGS] = TEMPLATE(INPUTARGS) Explain usage here
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% Outputs: 
% 	output1 - Description
% 	output2 - Description
% 
% Example: 
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015-2017
% Date: 15 August 2015
% Version: 0.1 (15 August 2015)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if numel(SYS.Room_Setup) > 1
        subSYS = SYS;
    end
    
    Broadband_Tools.Generate_Loudspeaker_Signals_batchfunc( subSYS );
    
end

