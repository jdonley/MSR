function [ d_framed ] = frame_data( d, OL, N )
%FRAME_DATA Frames data with overlap
% 
% Syntax:	[d_framed] = frame_data( d, OL, N ) 
% 
% Inputs: 
% 	d - Input data to frame
% 	OL - Overlap as a fraction (0.25 for 25% overlap)
% 	N - Frame length
% 
% Outputs: 
% 	d_framed - The framed input data
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
% Date: 16 August 2016
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find lengths of the overlapping and non-overlapping sections
no_OL_sz = floor(N*(1-OL));
OL_sz = N - no_OL_sz;

% Pad the signal file to allow the creation of the matrix
L = ceil( length(d) / no_OL_sz );
r = rem( length(d), no_OL_sz);
r(r==0)=no_OL_sz;
p = zeros( no_OL_sz - r , 1 );
d = [d; p]';

% Find the overlapping section and non-overlapping section
Half1  = reshape(d, [OL_sz L]).';
Half2 = [Half1(2:end, 1:no_OL_sz); zeros(1,OL_sz)];

% Concatenate the two sections
d_framed = [Half1 Half2 ];

% Remove last frame if it contains no extra information other than an overlap
if nnz(d_framed(end,:)) <= OL_sz
    d_framed = d_framed(1:(end-1),:);
    L = size(d_framed,1);
end

end

