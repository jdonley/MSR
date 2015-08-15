function [Q, R] = Gram_Schmidt2(A, w) %GRAM_SCHMIDT


    Q = A;
n_dimensions = size(A, 2);
R = zeros(n_dimensions);
R(1, 1) = norm(Q(:, 1));
Q(:, 1) = Q(:, 1) ./ R(1, 1);
for i = 2 : n_dimensions
    Qw = (Q(:, i - 1) .* w)' * Q(:, (i - 1) : end);
    R(i - 1, i : end) = Qw(2:end) / Qw(1);
    %% Surprisingly this loop beats the matrix multiply
    for j = i : n_dimensions
        Q(:, j) = Q(:, j) - Q(:, i - 1) * R(i - 1, j);
    end
    %% This multiply is slower than above
    %    Q(:, i : end) = ...
    %     Q(:, i : end) - ...
    %     Q(:, i - 1) * R(i - 1, i : end);
    R(i, i) = norm(Q(:,i));
    Q(:, i) = Q(:, i) ./ R(i, i);
end
    
end

