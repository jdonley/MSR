function [ Frames ] = PredictiveFraming( Frames, Buffer, Npredict, method )
%PREDICTIVEFRAMING This function accepts a framed signal and corresponding
%buffers to predict forward Npredict samples.
%
% Syntax:	[ Frames ] = PredictiveFraming( Frames, Buffer, Npredict, Npad, method )
%
% Inputs:
% 	Frames - Description
% 	Buffer - Description
% 	Npredict - Description
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
% Date: 16 August 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_p = Frames.'; % Frames_predicted
B = Buffer.';
[m,n1] = size(F_p);
[l,n2] = size(B);

if n1 ~= n2
    error('Number of frames mismatches with buffer.')
end

if strcmpi(method,'lpc')
    a = lpc(B,m-1);
elseif strcmpi(method,'burg')
    a = arburg(B,m-1);
elseif strcmpi(method,'cov')
    a = arcov(B,m/2);
end

A = -[zeros(n1,1) a(:,2:end)].';
y = zeros(Npredict,n1);
for f=1:n1-1
    [~, zf] = filter(A(:,f), 1, B(:,f));
    y(:,f) = filter([0 0], -a(f,:), zeros(1, Npredict), zf);
end
y_p = [B((end-m+Npredict+1):end,:); ...
    y;] ; %Predicted part
y_p(isnan(y_p))=0;

F_p = [F_p(m/2+1:end,:); ... %Latest part-frame known values
    y_p(1:m/2,:); ...
    y_p];

% F_p = [zeros(size(F_p,1),2), F_p(:,1:end-2)];

Frames = reshape(F_p,m,[]).';




end

