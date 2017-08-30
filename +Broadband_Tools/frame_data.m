function [ d_framed ] = frame_data( d, OL, N )
%Frames data with overlap
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
% See also: enframe

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016-2017
% Date: 28 August 2016
% Version: 0.2 (28 August 2017)
% Version: 0.1 (16 August 2016)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find lengths of the overlapping and non-overlapping sections
hopSz = double(int64(N*(1-OL)));    % Hop size
olSz = N - hopSz;                   % Overlap size

% Pad the signal file to allow the creation of the matrix
r = rem( length(d), hopSz);
r(r==0)=hopSz;
Lx = ceil(1/(1-OL));                % Number of extension frames
p = zeros( hopSz*Lx - r, 1 );           % Zero padding
d = [d; p].';                       % Extend length of input signal
L = length(d) / hopSz;              % Number of frames

% Find the overlapping section and non-overlapping section and concatenate
d_framed = reshape(d, [hopSz L]).';
OLinds = repmat((hopSz+1):N,L-Lx,1) ...
    + repmat((0:hopSz:hopSz*(L-Lx-1)).',1,olSz);
d_framed = [d_framed(1:L-Lx,:) d(OLinds)];
% d_framed = partFrames; 
% for seg = 1:
%     d_framed = [d_framed [partFrames((seg+1):end, 1:min(OL_sz,no_OL_sz)); zeros(seg,min(OL_sz,no_OL_sz))]];
% end

% Remove last frame if it contains no extra information other than an overlap
if nnz(d_framed(end,:)) <= olSz
    d_framed = d_framed(1:(end-1),:);
    L = size(d_framed,1);
end

end

