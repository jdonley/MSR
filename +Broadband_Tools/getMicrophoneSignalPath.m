function [ Path, Name, Ext, err, SubPath, spkr_sig_info_dir, Output_file_path_ext ] = getMicrophoneSignalPath( setup, signal_info, database_res, database_workingdir, method )
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
% Copyright: Jacob Donley 2016-2017
% Date: 17 August 2017
% Version: 0.2 (17 August 2017)
% Version: 0.1 (03 February 2016)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    method = 'new';
end

if nargin < 4
    database_workingdir = 'Z:\';
end

[ Path, Name, Ext, err, ...
    SubPath, spkr_sig_info_dir, Output_file_path_ext ] ...
    = Broadband_Tools.getSignalPath( ...
    setup, signal_info, database_res, database_workingdir, ...
    'Microphone_Signals', method );

end

