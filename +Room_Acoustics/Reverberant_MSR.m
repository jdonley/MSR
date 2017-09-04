function Reverberant_MSR(SYS)
% Summary of this function goes here
% 
% Syntax:	Reverberant_MSR(SYS)
% 
% Inputs: 
% 	SYS - Soundfield Reproduction system object
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
% Copyright: Jacob Donley 2016-2017
% Date: 14 June 2017
% Version: 0.1 (14 June 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose all;
delete(gcp('nocreate'));

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
    
    if numel(SYS.Room_Setup) > 1 ...
        && typ > 1
        % If there is more than one room and the first method has been 
        % completed then we choose the room set up for reproduction 
        % (transmission) and the associated loudspeaker setup
        subSYS = SYS;
        I_tx = strcmpi({subSYS.Room_Setup.SystemType},'transmit');
        subSYS.Main_Setup = SYS.Main_Setup(I_tx);
        subSYS.Room_Setup = SYS.Room_Setup(I_tx);
    end
    
    
    if isfield(SYS.signal_info, 'UseMeasuredATFs') && SYS.signal_info.UseMeasuredATFs, Room_Acoustics.useMeasuredATF(subSYS); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Room_Acoustics.Apply_RIRs.Reverberant_MSR_batchfunc( subSYS );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isfield(SYS.signal_info, 'UseMeasuredATFs') && SYS.signal_info.UseMeasuredATFs, Room_Acoustics.useSimulatedATF(subSYS); end % Revert back after use
    
end


