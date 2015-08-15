function [Q, R] = Gram_Schmidt_GPU(A, w, Q, R) %GRAM_SCHMIDT
%[m, n] = size(A);
%     Q  = complex(zeros(m, n,'like',A));
%     R  = complex(zeros(n, n,'like',A));
%     QQ = complex(zeros(m, n,'like',A));
%     %v  = zeros(n, 1,'like',A);
% 
%     for j = 1:n
%         v = A(:,j);
%         for i = 1:j-1
%             R(i,j) = (v.' * QQ(:,i));
%             v = v - R(i,j) * Q(:,i);
%         end
%         R(j,j) = norm(v);
%         Q(:,j) = v / R(j,j);
%         QQ(:,j) = (conj(Q(:,j)) .* w) ./ (w.' * (Q(:,j).*conj(Q(:,j))));
%     end

% A = gpuArray(A);    
% w = gpuArray(w);
    
    %Q = A;
     n_dimensions = size(A, 2);
%     R = zeros(n_dimensions ,'like',A);
    R(1, 1) = norm(Q(:, 1));
    Q(:, 1) = Q(:, 1) ./ R(1, 1);
    for i = 2 : n_dimensions
        Qw = (Q(:, i - 1) .* w)' * Q(:, (i - 1) : end);
        R(i - 1, i : end) = Qw(2:end) / Qw(1);
        
        Q(:, i : end) = ...
            Q(:, i : end) - ...
            Q(:, i - 1) * R(i - 1, i : end);
        R(i, i) = norm(Q(:,i));
        Q(:, i) = Q(:, i) ./ R(i, i);
    end
    
end

