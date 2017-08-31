function [Loudspeaker_Signals, Original_] = getPredictedLoudspeakerSignals( SYS )
%GETPREDICTEDLOUDSPEAKERSIGNALS Summary of this function goes here
% 
% Syntax:	[OUTPUTARGS] = GETPREDICTEDLOUDSPEAKERSIGNALS(INPUTARGS) Explain usage here
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
% Date: 31 August 2017 
% Version: 0.1 (31 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read Microphone Signals
MicSigPath = Broadband_Tools.getMicrophoneSignalPath( SYS );
MicSigFiles = Tools.keepFilesFromFolder( Tools.getAllFiles(MicSigPath), SYS.signal_info.speech_filepath);
MicSigFiles(~contains(MicSigFiles,'.mat'))=[];
MicSigFiles(~contains(MicSigFiles,SYS.signal_info.input_filename))=[];

MSF = load(MicSigFiles{:}); % MicSigFile (MSF)

%%


end
