function [DB, err, DB_path] = loadRIRDatabaseFromSetup( SYS_or_setup, room, database_workingdir, method )
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
% Date: 04 November 2015
% Version: 0.1 (04 November 2015)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SYS_type = 'Current_Systems.SR_System';

if nargin == 1
    if isa(SYS_or_setup,SYS_type)
        SYS = SYS_or_setup;        
        setup = SYS.Main_Setup;
        room = SYS.Room_Setup;
        database_workingdir = SYS.system_info.Drive;
    else
        error(['Single input argument must be of type: ' SYS_type]);
    end
else
    setup = SYS_or_setup;
end

%% 
if nargin < 4
    method = 'new2';
end
if nargin < 3
    database_workingdir = 'Z:\';
end

[DB_path, err] = Room_Acoustics.getRIRDatabasePath(setup,room,database_workingdir, method);

DB = load(DB_path);

end

