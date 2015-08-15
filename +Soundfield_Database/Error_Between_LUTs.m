clc;
clear;
close all;
tic; %Start timing this script
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))

%%
MSE_dBs = [];
Results_Path = '+Results\+3m_SpkrDia\';
SetupStyle = '65Spkrs_360DegArc';

N_f_low = 16;
N_w_low = 8;
N_increments = 5;
f = repmat(N_f_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1)) ,[N_increments 1]);
w = repmat(N_w_low,[N_increments N_increments]).*repmat(2.^(0:(N_increments-1))',[1 N_increments]);

i=1;j=0;
Input_file_path_NoLUT = [Results_Path '+' SetupStyle];

Input_file_paths = [repmat([Input_file_path_NoLUT '_LUT_'], [N_increments^2 1]) num2str(f(:)) repmat('f_',[N_increments^2 1]) num2str(w(:)) repmat('w\',[N_increments^2 1])];
Input_file_paths = strrep( num2cell(Input_file_paths,2), ' ', '');

Input_file_ext  = '.wav';
files = dir([Input_file_path_NoLUT '\*' Input_file_ext]);
%%
fprintf('\n====== Calculating LUT Comparison Statistics ======\n\n');
fprintf('\tCompletion: ');n=0;

for file = files'
    if strfind([Input_file_path_NoLUT '\' file.name(1:end-4) Input_file_ext],'Original') |  strfind([Input_file_path_NoLUT '\' file.name(1:end-4) Input_file_ext],'Bright')
        j=j+1;
    else
        x_NoLUT = audioread([Input_file_path_NoLUT '\' file.name(1:end-4) Input_file_ext]);
        x = [];
        for LUT = 1:(N_increments^2)
            x = audioread(cell2mat([Input_file_paths(LUT) file.name(1:end-4) Input_file_ext]));
            
            % MSE in dB
            %MSE_dBs(i, LUT) = pesq( ...
            %    [Input_file_path_NoLUT '\' file.name(1:end-4) Input_file_ext], ...
            %    cell2mat([Input_file_paths(LUT) file.name(1:end-4) Input_file_ext]));
            MSE_dBs(i, LUT) = mag2db( mean( (x - x_NoLUT).^2 ) );            
        end
        tElapsed = toc;
        ratio = (i+j)/length(files);
        tRem = (1-ratio) / ratio * tElapsed;
        tTot = tElapsed + tRem;
        fprintf(repmat('\b',1,n));
        n=fprintf('%.2f%% \n\tRemaining: %d mins %.0f secs \n\tTotal: %d mins %.0f secs\n', ratio * 100, floor(tRem/60), rem(tRem,60), floor(tTot/60), rem(tTot,60));
        i=i+1;
    end
end

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
ylim([min(MSE_dB_Means(:)-MSE_dB_95CI(:)) max(MSE_dB_Means(:)+MSE_dB_95CI(:))]);
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
save([Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__Statistics.mat'], ...
    'MSE_dB_Means', ...
    'MSE_dB_Stds', ...
    'MSE_dB_95CI');

saveas(h,[Results_Path '+LUT_Comparison_Statistics\' SetupStyle '__Statistics.fig']);
%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script
