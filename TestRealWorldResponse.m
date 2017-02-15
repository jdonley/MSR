clc;clear;close all;

fs = 16000;
fs_rec = 48000;
aliasFreq = (2*pi*24-1)*343 / (4*pi*0.9*pi);

% Simulated recording sweeps
recB = load('Z:\_RECOR~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_150HZ~1\SineSweepLin_80Hz-7500Hz_10s_Bright.mat');
recB = recB.Rec_Sigs_B';
recQ = load('Z:\_RECOR~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_150HZ~1\SineSweepLin_80Hz-7500Hz_10s_Quiet.mat');
recQ = recQ.Rec_Sigs_Q';

% Real world recording sweeps
recB_real = load('Z:\_RECOR~1\_PHYSI~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_1DB3B~1\SineSweepLin_80Hz-7500Hz_10s_Bright.mat');
recB_real = recB_real.Rec_Sigs_B;
recQ_real = load('Z:\_RECOR~1\_PHYSI~1\_CIRCL~1.3MP\_24GEN~1.261\_0BX_0~1.6QY\_VSRC_~1\_DATAB~1\_10X10~1\32REC_~1\_1DB3B~1\SineSweepLin_80Hz-7500Hz_10s_Quiet.mat');
recQ_real = recQ_real.Rec_Sigs_Q;


%% Resample real-world recordings
recB_real_16 = resample(double(recB_real),1,fs_rec/fs);
recQ_real_16 = resample(double(recQ_real),1,fs_rec/fs);

%% power norm signals
[~,scaleVal] = Broadband_Tools.power_norm(recB,recB_real_16,fs, [150, aliasFreq]);
recB_real_16 = recB_real_16 * scaleVal;
recQ_real_16 = recQ_real_16 * scaleVal;

%% Find power in frequency domain
Nfft = 2^(nextpow2(4000));
[pB_,f] = periodogram(recB,[],Nfft,fs);
pQ_ = periodogram(recQ,[],Nfft,fs);
pB = mean(pow2db(pB_),2);
pQ = mean(pow2db(pQ_),2);
[pB_real,f_real] = periodogram(recB_real_16,[],Nfft,fs);
pQ_real = periodogram(recQ_real_16,[],Nfft,fs);

f = f/1000; %to kHz
f_real = f_real/1000; %to kHz

%normalisation value
maxPow = max(pow2db([pB_real(:);pQ_real(:);pB_(:);pQ_(:)]));


%% Confidence Intervals
    pB_CI = Tools.confidence_intervals(pow2db(pB_'), 95);
    pQ_CI = Tools.confidence_intervals(pow2db(pQ_'), 95);

%% Generate confidence area plot data
pB_CI_area = ([pB-maxPow + pB_CI(:,1), ...
            pB_CI(:,2)-pB_CI(:,1) ]);
pQ_CI_area = ([pQ-maxPow + pQ_CI(:,1), ...
            pQ_CI(:,2)-pQ_CI(:,1) ]);

%% plot
brightCol = [0 0 1];
quietCol = [1 0 0];

figure(1);
plB_ = plot(f,pow2db(pB_)-maxPow,'color',1-(1-brightCol) * 0.2, 'linewidth',1); hold on;
plQ_ = plot(f,pow2db(pQ_)-maxPow,'color',1-(1-quietCol) * 0.2, 'linewidth',1); hold on;

arB = area(f, pB_CI_area,'LineStyle','none'); hold on;
arQ  = area(f, pQ_CI_area,'LineStyle','none');

plB = plot(f,pB-maxPow,'color',brightCol, 'linewidth',2); hold on;
plQ = plot(f,pQ-maxPow,'color',quietCol, 'linewidth',2); hold on;

plB_real = plot(f_real,pow2db(pB_real)-maxPow,'color',brightCol * 0.4, 'linewidth',2); hold on;
plQ_real = plot(f_real,pow2db(pQ_real)-maxPow,'color',quietCol * 0.4, 'linewidth',2); hold on;

% plot(f,pow2db(pB)-maxPow,'b', 'linewidth',2); hold on;
% plot(f,pow2db(pQ)-maxPow,'r', 'linewidth',2); hold on;



title('Frequency Responses');
ylabel('Power (dB)');
xlabel('Frequency (kHz)');
legend([plB_real,plQ_real,plB,plQ], ...
        {'Bright Zone Real-World (1-Mic)', ...
        'Quiet Zone Real-World (1-Mic)', ...
        'Bright Zone Simulated (32-Mics Mean)',...
        'Quiet Zone Simulated (32-Mics Mean)'}, ...
        'location', 'southeast');
set(gca,'XScale','log');
xlim([0.08 10]);
ylim([-80 0]);
grid minor

% for p = 1:length(plB_real)
% plB_real(p).Color = [plB_real(p).Color 0.2];
% plQ_real(p).Color = [plQ_real(p).Color 0.2];
% end

arB(1).FaceColor= 'none';
arQ(1).FaceColor= 'none';
arB(2).FaceColor= brightCol;
arQ(2).FaceColor= quietCol;
drawnow; pause(0.05);  % this is important for transparency!
arB(2).Face.ColorType = 'truecoloralpha';
arQ(2).Face.ColorType = 'truecoloralpha';
arB(2).Face.ColorData(4) = 0.4*255;
arQ(2).Face.ColorData(4) = 0.4*255;


hold off

%% plot 2
meanRange = [0.15 aliasFreq/1000];
f_meanrange = f>meanRange(1) & f<meanRange(2);
simMeanDiff = mean(pow2db(pB(f_meanrange)) - pow2db(pQ(f_meanrange)));
realMeanDiff = mean(pow2db(pB_real(f_meanrange)) - pow2db(pQ_real(f_meanrange)));

figure(2);
plot(f,pow2db(pB) - pow2db(pQ),'linewidth',2); hold on;
plot(f_real,pow2db(pB_real) - pow2db(pQ_real)); hold on;
plot(f_real,simMeanDiff+0*f_real, 'b','linewidth',2); hold on;
plot(f_real,realMeanDiff+0*f_real, 'r','linewidth',2); hold on;
hold off

text(7,simMeanDiff+4,['Mean = ' num2str(round(simMeanDiff,2)) 'dB'],'HorizontalAlignment','right');
text(7,realMeanDiff-4,['Mean = ' num2str(round(realMeanDiff,2)) 'dB'],'HorizontalAlignment','right');
text(10^mean(log10(meanRange)),-10, ...
    ['Frequency Range for Mean: ' ...
    num2str(round(meanRange(1),2)) 'kHz to ' ...
    num2str(round(meanRange(2),2)) 'kHz'], ...
    'HorizontalAlignment','center');

title('Level Differences (Suppression)');
ylabel('Power (dB)');
xlabel('Frequency (kHz)');
legend({'Suppression Simulated (32-Mics)', ...
        'Suppression Real-World (1-Mic)', ...
        'Mean Simulated Suppression (32-Mics)', ...
        'Mean Real-World Suppression (1-Mic)'}, ...
        'location', 'northeast');
hold on;plot(f_real,0*f_real, 'k','linewidth',2);hold off; %plot zero
hold on;plot(meanRange(1)*ones(2,1),[-20,80], '--k','linewidth',1);hold off; %plot mean range
hold on;plot(meanRange(2)*ones(2,1),[-20,80], '--k','linewidth',1);hold off; %plot mean range
set(gca,'XScale','log');
xlim([0.08 8]);
ylim([-20 80]);
grid minor

