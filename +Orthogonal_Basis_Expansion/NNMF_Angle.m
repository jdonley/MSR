function [W, H] = NNMF_Angle(A, w) %GRAM_SCHMIDT
[m, n] = size(A);
    W_  = complex(zeros(m, n-1));
    H_  = complex(zeros(n-1, n));
    W  = complex(zeros(m, n-1));
    H  = complex(zeros(n-1, n));
    
R = real(A);
R = R - min(R(:));
[W_, H_] = nnmf(R, n-1,'algorithm','mult');%,'replicate',5);
W = W + W_;
H = H + H_;

I = imag(A);
I = I - min(I(:));
[W_, H_] = nnmf(I, n-1,'algorithm','mult');%,'replicate',5);
W = W + W_ * 1i;
H = H + H_ * 1i;



%W = W - pi;
%W = exp( 1i * W );

W(:,n) = zeros(m,1);
H(n,:) = zeros(1,n);
end

