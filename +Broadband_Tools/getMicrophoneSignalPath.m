function [ Path, Name, err ] = getMicrophoneSignalPath( SYS )
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
% Copyright: Jacob Donley 2017
% Date: 17 August 2017
% Version: 0.1 (17 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Path, err ] = Results.getRecordingsPath( SYS );
Name     = [SYS.signal_info.input_filename SYS.system_info.sc ];

end

