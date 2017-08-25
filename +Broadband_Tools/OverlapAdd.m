function [ TimeSignal ] = OverlapAdd( FramedSig, Overlap )
% Performs overlap-add on the framed time signal
% 
% Syntax:	[ TimeSignal ] = OverlapAdd( FramedSig, Overlap )
% 
% Inputs: 
% 	FramedSig - The framed signal as a 2D matrix (Frames x Time)
% 	Overlap - The amount of overlap as a decimal number (0.5 = 50% overlap)
% 
% Outputs: 
% 	TimeSignal - The one dimensional time domain signal
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
% Date: 15 August 2015
% Version: 0.1 (15 August 2015)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

one_minus_overlap_size = floor(size(FramedSig,2)*(1-Overlap));
overlap_size = size(FramedSig,2) - one_minus_overlap_size;

first_few = FramedSig(1, 1:overlap_size);
last_few  = FramedSig(end, one_minus_overlap_size+1:end);

overlapping_1 = FramedSig(1:end-1, one_minus_overlap_size+1:end);
overlapping_2 = FramedSig(2:end  , 1:overlap_size  );

added_part = (overlapping_1 + overlapping_2)';

non_overlapping = FramedSig(:, overlap_size+1:one_minus_overlap_size);

if ~isempty(non_overlapping)
    TimeSignal = [non_overlapping(1:end-1,:) added_part']';
    TimeSignal = TimeSignal(:);
    TimeSignal = [first_few'; TimeSignal; non_overlapping(end,:)'; last_few'];
else
    TimeSignal = [first_few'; added_part(:); last_few'];
end

end

