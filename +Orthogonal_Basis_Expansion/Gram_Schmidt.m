function [Q, R] = Gram_Schmidt(A, w)
% Gram-Schmidt orthogonalisation with inner-product weighting
% 
% Syntax:	[Q, R] = Gram_Schmidt(A, w)
% 
% Inputs: 
% 	A - Set of vectors to orthonormalise
% 	w - Inner-product weights
% 
% Outputs: 
% 	Q - Orthogonal matrix
% 	R - Triangluar matrix
% 
% Example: 
%     A = complex( rand(360000,100), rand(360000,100));
%     w = rand(360000, 1);
%     [Q, R] = Gram_Schmidt(A, w);
% 
% See also: norm, conj, complex

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015-2017
% Date: 15 August 2015
% Version: 0.1 (15 August 2015)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n] = size(A);
    Q  = complex(zeros(m, n));
    R  = complex(zeros(n, n));
    QQ = complex(zeros(m, n));

    for j = 1:n
        v = A(:,j);
        for i = 1:j-1
            R(i,j) = (v.' * QQ(:,i));
            v = v - R(i,j) * Q(:,i);
        end
        R(j,j) = norm(v);
        Q(:,j) = v / R(j,j);
        QQ(:,j) = (conj(Q(:,j)) .* w) ./ (w.' * (Q(:,j).*conj(Q(:,j))));
    end

end

