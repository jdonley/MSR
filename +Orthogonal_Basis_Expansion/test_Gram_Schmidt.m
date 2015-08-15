clc;
clear;
close all;

N = 600
samples = (ceil(N^(1/2)) + 1)^2
width = ceil(N^(1/2)) + 1
centre = ceil(width/2)

 A=complex( rand(samples,N), rand(samples,N));
 w=rand(samples,1);
% AA = gpuArray(A);
% ww = gpuArray(w);
% Q2 = AA;
% R2 = zeros(size(AA, 2) ,'like',AA);

% A=complex( rand(100,10), rand(100,10));
% w=rand(100,1);
tic;
[Q, R] = Orthogonal_Basis_Expansion.Gram_Schmidt(A, w);
toc;
tic;
[Q2, R2] = Orthogonal_Basis_Expansion.Gram_Schmidt_mex(A, w);
toc;
% tic;
% [Q2, R2] = Orthogonal_Basis_Expansion.Gram_Schmidt_GPU(AA, ww, Q2, R2);
% toc

norm(Q-Q2)
norm(R-R2)

%rcond(R)
%rcond(R2)