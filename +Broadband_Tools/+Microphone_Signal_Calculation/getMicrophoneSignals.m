function [ Mic_Signals, Mic_Signals_Anechoic, Original ] = getMicrophoneSignals( Input_Signal, SYS )
% Summary of this function goes here
% 
% Syntax:	[ Mic_Signals, Original_ ] = getMicrophoneSignals( Input_Signal, SYS )
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

DB = Room_Acoustics.loadRIRDatabaseFromSetup(SYS);

Mic_Signals = Tools.fconv( ...
    repmat(Input_Signal,1,SYS.Room_Setup.NoReceivers), ...
    DB.RIRs.Bright_RIRs.');

if isfield(DB.RIRs,'Bright_RIRs_Anechoic') 
    Mic_Signals_Anechoic = Tools.fconv( ...
        repmat(Input_Signal,1,SYS.Room_Setup.NoReceivers), ...
        DB.RIRs.Bright_RIRs_Anechoic.');
else
    Mic_Signals_Anechoic = [];
end

Original = Input_Signal;

end

