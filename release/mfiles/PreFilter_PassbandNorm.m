function [ y_Leveled ] = PreFilter_PassbandNorm( y_orig, y_diff, setup, signal_info )
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
% Date: 7 December 2016 
% Revision: 0.3
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = numel(y_orig(:));
if ~mod(N,2) % If N is even
    n = 0.5;
else % If N is odd
    n = (N-1)/(2*N);
end

f_cutoff = min( [...
    Broadband_Tools.getAliasingFrequency(setup) * signal_info.c / (2*pi), ...
    n * signal_info.Fs ] );

y_normed = y_orig(:) ./ rms(y_orig(:)) * db2mag( signal_info.clipLevelAdjust ); 
y_Leveled = Broadband_Tools.power_norm( y_normed, y_diff(:), signal_info.Fs, [signal_info.f_low, f_cutoff]);

end

