clear;clc;
fs = 16000;

prepath = 'Z:\+Recordings\+LINEarray_1.3mPerpDist_180degCentre\';

load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_-0.6By_0Qx_0.6Qy\+VSrc_ps_180deg_1.3m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_0dB_10000weight__method_ZoneWeightMaskerAliasCtrl\Masker_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_-0.6By_0Qx_0.6Qy\+VSrc_ps_180deg_1.3m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_0dB_10000weight__method_ZoneWeightMaskerAliasCtrl\Masker_Quiet.mat'])
Rec_Sigs_B_m = Rec_Sigs_B;
Rec_Sigs_Q_m = Rec_Sigs_Q;
load([ prepath '+1ParametricSpkrs_0mLen\+0Bx_-0.6By_0Qx_0.6Qy\+VSrc_ps_180deg_1.3m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_0dB_1weight__method_ParametricMaskerAliasCtrlHPF\Masker_Bright.mat'])
load([ prepath '+1ParametricSpkrs_0mLen\+0Bx_-0.6By_0Qx_0.6Qy\+VSrc_ps_180deg_1.3m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_0dB_1weight__method_ParametricMaskerAliasCtrlHPF\Masker_Quiet.mat'])
Rec_Sigs_B_p = Rec_Sigs_B;
Rec_Sigs_Q_p = Rec_Sigs_Q;


load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female_SA2_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female_SA2_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male_SA2_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male_SA2_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female_SX126_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female_SX126_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male_SI1669_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male_SI1669_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female_SX306_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female_SX306_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female1_SI1544_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Female1_SI1544_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male_SX49_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male_SX49_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male1_SX5_Bright.mat'])
load([ prepath '+24Genelec8010ASpkrs_2.806mLen\+0Bx_0.6By_0Qx_-0.6Qy\+VSrc_pw_0deg_0m\+Database_512f_32w\+10x10x10Dim_1Ab\32Rec_5x5x5Ctr\+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\Male1_SX5_Quiet.mat'])
Rec_Sigs_B_sp = Rec_Sigs_B;
Rec_Sigs_Q_sp = Rec_Sigs_Q;
pxx_b_sp(:,end+1:end+32) = periodogram(Rec_Sigs_B_sp',[],1024,fs);
pxx_q_sp(:,end+1:end+32) = periodogram(Rec_Sigs_Q_sp',[],1024,fs);

%%
gain = -15;%db

b_mask = Rec_Sigs_B_m(:,1:length(Rec_Sigs_B_p)) + Rec_Sigs_B_p;
q_mask = Rec_Sigs_Q_m(:,1:length(Rec_Sigs_Q_p)) + Rec_Sigs_Q_p;

 b_mask = Rec_Sigs_B_m(:,1:length(Rec_Sigs_B_p));
 q_mask = Rec_Sigs_Q_m(:,1:length(Rec_Sigs_Q_p));
% b_mask = Rec_Sigs_B_p;
% q_mask = Rec_Sigs_Q_p;


pxx_b_mask = periodogram(b_mask',[],1024,fs) .* db2mag(gain);
pxx_q_mask = periodogram(q_mask',[],1024,fs) .* db2mag(gain);


f = linspace(0,fs/2,length(pxx_b_mask));
pxx_b_mask_smooth = smooth(mag2db(mean(pxx_b_mask,2)),20);
pxx_q_mask_smooth = smooth(mag2db(mean(pxx_q_mask,2)),20);
pxx_b_sp_smooth = smooth(mag2db(mean(pxx_b_sp,2)),20);
pxx_q_sp_smooth = smooth(mag2db(mean(pxx_q_sp,2)),20);

%%

figure(1); hold off;
plot(f, (pxx_b_mask_smooth) ,'b.'); hold on
plot(f, (pxx_q_mask_smooth) ,'r.');
plot(f, (pxx_b_sp_smooth) ,'b');
plot(f, (pxx_q_sp_smooth) ,'r');

legend({ ...
        'Masker in Bright';
        'Masker in Quiet';
        'Speech in Bright';
        'Speech in Quiet'});

set(gca,'XScale','log');
xlim([100 10000]);
ylim([-280 -100]);
grid minor

figure(2);hold off
plot(f, (pxx_b_sp_smooth - pxx_b_mask_smooth) ,'b'); hold on
plot(f, (pxx_q_mask_smooth - pxx_q_sp_smooth) ,'r');
avg1 = repmat(mean(pxx_b_sp_smooth - pxx_b_mask_smooth),size(pxx_b_sp_smooth));
avg2 = repmat(mean(pxx_q_mask_smooth - pxx_q_sp_smooth),size(pxx_b_sp_smooth));
plot(f, avg1 ,'b');
plot(f, avg2 ,'r');

plot(f, (((pxx_b_sp_smooth - pxx_b_mask_smooth) - avg1 <0) & ((pxx_q_mask_smooth - pxx_q_sp_smooth) - avg2 > 0)) .* -10 ,'k');
plot(f, (((pxx_b_sp_smooth - pxx_b_mask_smooth) - avg1 >0) & ((pxx_q_mask_smooth - pxx_q_sp_smooth) - avg2 < 0)) .* 10 ,'k');


set(gca,'XScale','log');
xlim([100 10000]);
grid minor