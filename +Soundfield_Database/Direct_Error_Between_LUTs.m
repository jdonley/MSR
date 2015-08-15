clc;
clear;
close all;
tic; %Start timing this script
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%%
MSE_dBs = [];
Results_Path = '+Results\+3m_SpkrDia\';
Database_Path = '+Soundfield_Database\+3m_SpkrDia\';
SetupStyle = '65Spkrs_360DegArc';
interpMethod = 'bilinear';

pw_angle = 15;
Fs = 16000;
Nfft = 1024;
Frequencies_ = linspace(0, Fs/2, Nfft/2+1);
Frequencies_ = Frequencies_(1:end-1);
Weights_ = [0     logspace(log10( 1e-2), log10(  1e4), 256 - 1) ];

N_f_low = 16;
N_w_low = 8;
N_increments = 5;
f = repmat(N_f_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1)) ,[N_increments 1]);
w = repmat(N_w_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1))',[1 N_increments]);
f_ = f(:);
w_ = w(:);

%%
fprintf('\n====== Calculating LUT Comparison Statistics ======\n\n');
fprintf('\tCompletion: ');n=0;

for LUT = 1:(N_increments^2)
    load([Database_Path '+' SetupStyle '\LUT_Weight_vs_Frequency_' num2str(pw_angle) 'deg_' num2str(f_(LUT)) 'f_' num2str(w_(LUT)) 'w' '.mat']);
    LUT_Bright = Bright_Sample__Weight_Vs_Frequency;%(db2mag(Contrast__Weight_Vs_Frequency).^0.5 .* db2mag(Quiet_SPL__Weight_Vs_Frequency-94));
    LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;%db2mag(Quiet_SPL__Weight_Vs_Frequency-94);
    
%     for weight=1:length(Weights_)        
%         Bright_Zone_sample2(:,weight,LUT) = permute( Tools.interpVal_2D(LUT_Bright, Frequencies, Weights, Frequencies_, Weights_(weight)), [2 1]);
%         Quiet_Zone_sample2(:,weight,LUT)  = permute( Tools.interpVal_2D(LUT_Quiet, Frequencies, Weights, Frequencies_, Weights_(weight)), [2 1]);
%     end
    Bright_Zone_sample(:,:,LUT)= imresize(LUT_Bright,[w_(N_increments^2)*2 f_(N_increments^2)*2], interpMethod);
    Quiet_Zone_sample(:,:,LUT) = imresize(LUT_Quiet,[w_(N_increments^2)*2 f_(N_increments^2)*2], interpMethod);
%     if LUT == 1
%         mean_b = Bright_Zone_sample(:,:,LUT);
%         mean_b(isnan(mean_b)) = 0;
%         mean_q = Quiet_Zone_sample(:,:,LUT);
%         mean_q(isnan(mean_q)) = 0;
%         Bright_Zone_sample(:,:,LUT) = ones(size(Bright_Zone_sample(:,:,LUT))) * mean(mean_b(:));
%         Quiet_Zone_sample(:,:,LUT)  = ones(size(Quiet_Zone_sample (:,:,LUT))) * mean(mean_q(:));
%     end
    
    tElapsed = toc;
    ratio = LUT / N_increments^2;
    tRem = (1-ratio) / ratio * tElapsed;
    tTot = tElapsed + tRem;
    fprintf(repmat('\b',1,n));
    n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
    
end

load([Database_Path '+' SetupStyle '\LUT_Weight_vs_Frequency_' num2str(pw_angle) 'deg_' num2str(f_(N_increments^2)*2) 'f_' num2str(w_(N_increments^2)*2) 'w' '.mat']);
LUT_Bright = Bright_Sample__Weight_Vs_Frequency;%(db2mag(Contrast__Weight_Vs_Frequency).^0.5 .* db2mag(Quiet_SPL__Weight_Vs_Frequency-94));
LUT_Quiet  = Quiet_Sample__Weight_Vs_Frequency;%db2mag(Quiet_SPL__Weight_Vs_Frequency-94);

% for weight=1:length(Weights_)
%     Bright_Zone_ref(:,weight) = permute( Tools.interpVal_2D(LUT_Bright, Frequencies, Weights, Frequencies_, Weights_(weight)), [2 1]);
%     Quiet_Zone_ref(:,weight)  = permute( Tools.interpVal_2D(LUT_Quiet, Frequencies, Weights, Frequencies_, Weights_(weight)), [2 1]);
% end
Bright_Zone_ref = imresize(LUT_Bright,[w_(N_increments^2)*2 f_(N_increments^2)*2], interpMethod);
Quiet_Zone_ref  = imresize(LUT_Quiet,[w_(N_increments^2)*2 f_(N_increments^2)*2], interpMethod);

for LUT = 1:(N_increments^2)
    diff_b = Bright_Zone_ref' - Bright_Zone_sample(:,:,LUT)';
    diff_b(isnan(diff_b))=0;
    
    diff_q = Quiet_Zone_ref' - Quiet_Zone_sample(:,:,LUT)';
    diff_q(isnan(diff_q))=0;
    
   MSE_dBs(:,LUT) =  [diff_b(:); diff_q(:)] .^ 2;   
end
%MSE_dBs(:,N_increments^2) =  zeros(size(MSE_dBs,1),1);

MSE_dBs = mag2db(MSE_dBs);
MSE_dBs(MSE_dBs == -Inf) = 0;

%%
MSE_dBs_temp = reshape( MSE_dBs(MSE_dBs ~= -Inf), [], N_increments^2);

MSE_dB_Means = mean(MSE_dBs_temp, 1);
MSE_dB_Stds = std(MSE_dBs_temp, 0, 1); %Sample standard deviation (N-1) caused by flag=0

p = 0.95;
z = sqrt(2)*erfcinv((1-p));
MSE_dB_95CI = z * MSE_dB_Stds / sqrt(size(MSE_dBs_temp,1));

MSE_dB_Means = reshape(MSE_dB_Means, N_increments, [])';
MSE_dB_Stds = reshape(MSE_dB_Stds, N_increments, [])';
MSE_dB_95CI = reshape(MSE_dB_95CI, N_increments, [])';

%%
h = bar(MSE_dB_Means);
for s = 1:N_increments
    set(h(s), 'EdgeColor', get(h(s),'FaceColor'));
    set(h(s), 'FaceColor', 'white');
    set(h(s), 'LineWidth', 3);
end
ylim([min(MSE_dB_Means(:)-MSE_dB_95CI(:)) -70]);
xlabel('LUT No of Frequencies');
ylabel('Mean Squared Error (dB)');
set(gca,'YDir','reverse');
set(h,'BarWidth',1);    % The bars will now touch each other
set(gca,'YGrid','on')
set(gca,'GridLineStyle','-');
set(gca, 'XTickLabel', f(1,:));
l = legend(num2str(w(:,1)), 'Location','northeastoutside'); t = get(l,'title');
set(t, 'string', 'LUT No of Weights');
hold on;
numgroups = size(MSE_dB_Means, 1);
numbars = size(MSE_dB_Means, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    h = errorbar(x, MSE_dB_Means(:,i), MSE_dB_95CI(:,i), 'k', 'linestyle', 'none', 'LineWidth', 2);
end
hold off;

%%
if ~exist([Results_Path '+LUT_Comparison_Statistics\'],'dir'); mkdir([Results_Path '+LUT_Comparison_Statistics\']); end
%save([Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__Statistics.mat'], ...
save([Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__LUT2LUT_Statistics.mat'], ...
    'MSE_dB_Means', ...
    'MSE_dB_Stds', ...
    'MSE_dB_95CI');

%saveas(h,[Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__Statistics.fig']);
saveas(h,[Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__LUT2LUT_Statistics.fig']);

path=[Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__LUT2LUT_Statistics.xlsx'];
xlswrite(path, MSE_dB_Means,1,'C3');
xlswrite(path, MSE_dB_Stds,2,'C3');
xlswrite(path, MSE_dB_95CI,3,'C3');
for sheet = 1:3    
xlswrite(path, {'Weights'},sheet,'C1');
xlswrite(path, w(:,1)',sheet,'C2');
xlswrite(path, {'Frequencies'},sheet,'A3');
xlswrite(path, f(1,:)',sheet,'B3');
end
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
