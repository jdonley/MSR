clear
Fs = 16000;
Nfft = 1024;
Frequencies_ = linspace(0, Fs/2, Nfft/2+1);
Frequencies_ = Frequencies_(1:end-1);
Weights_ = [0     logspace(log10( 1e-2), log10(  1e4), 128 - 1) ];

N_f_low = 16;
N_w_low = 8;
N_increments = 5;
f = repmat(N_f_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1)) ,[N_increments 1]);
w = repmat(N_w_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1))',[1 N_increments]);
f = f(:);
w = w(:);

for i = 1:length(f)
    load(['+Soundfield_Database\+' num2str(1.5*2) 'm_SpkrDia\+' num2str(65) 'Spkrs_' num2str(360) 'DegArc\LUT_Weight_vs_Frequency_' num2str(15) 'deg_' num2str(f(i)) 'f_' num2str(w(i)) 'w' '.mat']);
    LUT_Bright = Bright_Sample__Weight_Vs_Frequency;%(db2mag(Contrast__Weight_Vs_Frequency).^0.5 .* db2mag(Quiet_SPL__Weight_Vs_Frequency-94));
    LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;%db2mag(Quiet_SPL__Weight_Vs_Frequency-94);
    
    for j=1:length(Weights_)
        Bright_Zone_sample(:,j,i) = permute( Tools.interpVal_2D(LUT_Bright, Frequencies, Weights, Frequencies_, Weights_(j)), [2 1]);
        Quiet_Zone_sample(:,j,i)  = permute( Tools.interpVal_2D(LUT_Quiet, Frequencies, Weights, Frequencies_, Weights_(j)), [2 1]);
    end
    %figure(i);surf(Frequencies_, Weights_, mag2db(Quiet_Zone_sample(:,:,i)'),'LineStyle','none'); view(2);set(gca,'XScale','log');set(gca,'YScale','log');

end

for  i = 1:(length(f)-1)
    diff = Bright_Zone_sample(:,:,end)' - Bright_Zone_sample(:,:,i)';
    diff(isnan(diff))=0;
   Difference(i) =  mean(mean((diff .^ 2)));
   diff = Quiet_Zone_sample(:,:,end)' - Quiet_Zone_sample(:,:,i)';
    diff(isnan(diff))=0;
   Difference(i) =  mean([Difference(i) mean(mean((diff .^ 2)))]);
end
mag2db(Difference)'