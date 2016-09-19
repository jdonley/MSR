function [ Frames ] = PredictivePadding( Frames, Buffer, Npredict, Npad, method )
%PREDICTIVEPADDING Summary of this function goes here
%
% Syntax:	[OUTPUTARGS] = TEMPLATE(INPUTARGS) Explain usage here
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
% Date: 16 August 2016
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_p = Frames.'; % predicted padding
B = Buffer.';
[m,n1] = size(F_p);
[l,n2] = size(B);

if n1 ~= n2
    error('Frame number mismatch with buffer.')
end

if strcmpi(method,'lpc')
    a = lpc(B,m-1);
elseif strcmpi(method,'burg')
    a = arburg(B,m-1);
end

A = -[zeros(n1,1) a(:,2:end)].';
y = zeros(Npredict,n1);
for f=1:n1
    [~, zf] = filter(A(:,f), 1, B(:,f));
    y(:,f) = filter([0 0], -a(f,:), zeros(1, Npredict), zf);
end
F_p = [F_p; ... %Known part
         y; ... %Predicted part
    zeros(Npad - Npredict, size(F_p,2))]; %Zero padding

% Frames = [zeros(size(F_p,2),Npad/2), ...
%     F_p.', ...
%     zeros(size(F_p,2),Npad/2 - Npredict)];

Frames = circshift(F_p,Npad/2,1).';


end

