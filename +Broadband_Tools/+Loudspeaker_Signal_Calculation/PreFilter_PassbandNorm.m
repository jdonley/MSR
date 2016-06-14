function [ y_Leveled ] = PreFilter_PassbandNorm( y_orig, y_diff, setup, signal_info, levelAdjust )
%PREFILTER_PASSBANDNORM Adjust the signal to account for the aliasing caused by a limited number of loudspeakers
% 
% Syntax:	[OUTPUTARGS] = ARBITRARYOCTAVEFILT(INPUTARGS) Explain usage here
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
if nargin < 5
    levelAdjust = -35; % Hope that -35dB RMS level is low enough to avoid clipping upon saving
end

f_cutoff = Broadband_Tools.getAliasingFrequency(setup) * signal_info.c / (2*pi);

y_normed = y_orig(:) ./ rms(y_orig(:)) * db2mag(levelAdjust); 
y_Leveled = Broadband_Tools.power_norm( y_normed, y_diff(:), signal_info.Fs, [signal_info.f_low, f_cutoff]);

end

