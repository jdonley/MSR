c = 343;
fs = 16000;
cTs = c/fs;

LL=[4 9 3];
LL=[8 10 3];
LL=[9 14 3];
%LL = [10 15 12.5].*0.3048;
nsamples=8000;

L(1) = LL(1)/cTs;
L(2) = LL(2)/cTs;
L(3) = LL(3)/cTs;
n1 =  ceil(nsamples/(2*L(1)));
n2 =  ceil(nsamples/(2*L(2)));
n3 =  ceil(nsamples/(2*L(3)));
paths = 8*(2*n1+1)*(2*n2+1)*(2*n3+1)
maxorder = abs(2*n1 - 0) + abs(2*n2 - 0) + abs(2*n3 - 0)


% c = 340;                    % Sound velocity (m/s)
% fs = 16000;                 % Sample frequency (samples/s)
% r = [1.5 2 3];              % Receiver position [x y z] (m)
% s = [2.5 2 3];              % Source position [x y z] (m)
% L = [4 9 3];                % Room dimensions [x y z] (m)
% beta = 0.3;                 % Reverberation time (s)
% n = 8000;                   % Number of samples
% mtype = 'omnidirectional';
% order = -1;
% 
% h = rir_generator(c, fs, r, s, L, beta, n, mtype, order);
% h2= rir_generator(c, fs, r, s, L, beta, n, mtype, 74);
% norm(h-h2)