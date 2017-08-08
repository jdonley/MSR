function useSimulatedATF(SYS)
%USESIMULATEDATF Summary of this function goes here
% 
% Syntax:	[OUTPUTARGS] = USESIMULATEDATF(INPUTARGS) Explain usage here
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
% Copyright: Jacob Donley 2017
% Date: 26 April 2017 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    SYS = Current_Systems.loadCurrentSRsystem;
end

%% Find correct setups
ArrType = SYS.system_info.CurrentSpeakerArrayType;
% isCircArray = contains(ArrType, 'circ');
isLineArray = contains(ArrType, 'line');
I = cellfun(@(x) contains(x,ArrType) ,{SYS.Main_Setup.Speaker_Array_Type},'un',0);
S = SYS.Main_Setup([I{:}]);

if isempty(S), return; end

%% Load current RIRs
[DB,~,RIR_FilePath] = Room_Acoustics.loadRIRDatabaseFromSetup( S(1), SYS.Room_Setup, SYS.system_info.Drive );
if ~(isfield(DB,'isRecordedRIR') && ~DB.isRecordedRIR)
    
    % Get path strings
    fpath = fileparts(RIR_FilePath);
    RIRfiles = Tools.getAllFiles(fpath);
    DBs = cellfun(@load , RIRfiles, 'un', 0);
    
    hasIsRecField = cellfun(@(d) isfield(d,'isRecordedRIR'), DBs,'un',0);
    hasIsRecField = [hasIsRecField{:}];
    DBsWithIsRecFld = [DBs{hasIsRecField}];
    SimRIRs = hasIsRecField;
    SimRIRs(hasIsRecField) = ~[DBsWithIsRecFld.isRecordedRIR];
    
    if any(SimRIRs)
        DB = DBs{SimRIRs};
        DB = DB(1); % Only use one if there are multiple
    else
        error('No simulated RIR databases exist');
    end
end
RIRs = DB.RIRs;
isRecordedRIR = DB.isRecordedRIR;

%% Re-save simulated RIRs
save(RIR_FilePath,'RIRs','isRecordedRIR');

end
