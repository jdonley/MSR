p=parpool;

%%
fs = 16000;
spkrs = 295;
recs = 32;
sig_len = fs*2.5;
rir_len = 8000;

Speaker_Signals = rand( spkrs,sig_len );
RIRs = rand( recs, rir_len, spkrs );

%%
Receiver_Signals1 = zeros( spkrs, recs, sig_len + rir_len - 1 );
Receiver_Signals2 = zeros( spkrs, recs, sig_len + rir_len - 1 );
%Receiver_Signals3 = gpuArray.zeros( spkrs, recs, sig_len + rir_len - 1 );

%%
tic;
Receiver_Signals1 = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs(Speaker_Signals,RIRs, 'FFT_slow');
toc

tic;
Receiver_Signals2 = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs(Speaker_Signals,RIRs, 'FFT');
toc

% tic;
% Receiver_Signals3 = Room_Acoustics.Apply_RIRs.Convolve_SpkrSigs_and_RIRs(Speaker_Signals,RIRs, 'GPU');
% toc

%%
norm(Receiver_Signals1(:)-Receiver_Signals2(:))
%norm(Receiver_Signals2(:)-Receiver_Signals3(:))
%norm(Receiver_Signals1(:)-Receiver_Signals3(:))

%%
delete(p);