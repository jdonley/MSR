function [ y_PreFilt ] = PreFilter_AliasCtrl( y, setup, signal_info, cheby_order, ripple )
%ALIASCTRL_PREFILTER Adjust the signal to account for the aliasing caused by a limited number of loudspeakers
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
% Date: 10 December 2016
% Revision: 0.3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    cheby_order=9;
    ripple = 1;
end
f_cutoff = Broadband_Tools.getAliasingFrequency(setup) * signal_info.c / (2*pi);

if f_cutoff < (signal_info.Fs/2)
    % Create filter
    [b,a] = cheby1(cheby_order, ripple, f_cutoff./(signal_info.Fs/2), 'low');
    
    if isstable(b,a) % Filter the signal if the filter is stable
        % Apply low pass filter to signal
        y_PreFilt = filter(b,a,y(:));
    else % If the filter is not stable return the original signal
        % (filtering is pointless for fringe cases)
        y_PreFilt = y(:);
    end
else % If the cutoff is not within the interval (0,1) return the original signal
    y_PreFilt = y(:);
end

end

