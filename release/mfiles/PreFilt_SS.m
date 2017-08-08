function [y_SS, spect, frqs] = PreFilt_SS( y, signal_info )
%PREFILT_SS Summary of this function goes here
% 
% Syntax:	[OUTPUTARGS] = PREFILT_SS(INPUTARGS) Explain usage here
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
% Copyright: Jacob Donley 2016
% Date: 06 June 2016 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[spect,frqs]=Tools.LTASS( signal_info.speech_filepath, signal_info.Nfft );

y_SS = Tools.ArbitraryOctaveFilt(y, spect, frqs, signal_info.Nfft, signal_info.Fs, signal_info.OctaveBandSpace);

end
