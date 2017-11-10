function [ sig_matched, adjustVal ] = power_norm( sig_orig, sig_diff, Fs, band)
% Summary of this function goes here
% 
% Syntax:	[ sig_matched, adjustVal ] = ...
%               power_norm( sig_orig, sig_diff, Fs, band)
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
% Date: 22 February 2016
% Version: 0.1 (22 February 2016)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(band) == 1
    % Calculate octave band
    band = band * [2^-.5 2^.5];
end

bpow_in   = mean(bandpower( sig_orig, Fs, band));
bpow_diff = mean(bandpower( sig_diff, Fs, band));
adjustVal = mean( db2mag( pow2db( bpow_in ./ bpow_diff ) ) );
sig_matched = sig_diff * adjustVal;

end

