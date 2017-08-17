function Generate_Microphone_Signals(SYS)
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
% Date: 17 August 2017
% Version: 0.1 (17 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load System
if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end

%%

for r = 1:numel(SYS.Room_Setup)
    
    if strcmpi(SYS.Room_Setup(r).SystemType, 'receive')
        subSYS = SYS;
        subSYS.Main_Setup = SYS.Main_Setup(r);
        subSYS.Room_Setup = SYS.Room_Setup(r);
        
        Broadband_Tools.Generate_Microphone_Signals_batchfunc( subSYS );
    end
end

