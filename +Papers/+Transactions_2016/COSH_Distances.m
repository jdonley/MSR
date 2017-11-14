clear;
%%
SYS = Current_Systems.IEEETransactions_System_D;
SYS.Main_Setup = SYS.Main_Setup(2);
SYS.Masker_Setup = SYS.Masker_Setup(2);
res_path = Results.getRecordingsPath(SYS);

files = Tools.getAllFiles( res_path );
files = Tools.keepFilesFromFolder(files,SYS.signal_info.speech_filepath);
filesQ = files(find(~cellfun(@isempty,strfind( files, 'Quiet'))));
filesB = files(find(~cellfun(@isempty,strfind( files, 'Bright'))));
speech = [];speechB = [];
for f = 1:numel(filesQ)
    S = load(filesQ{f});
    speech = [speech; S.Rec_Sigs_Q.'];
end
for f = 1:numel(filesB)
    S = load(filesB{f});
    speechB = [speechB; S.Rec_Sigs_B.'];
end
win_=rectwin(SYS.analysis_info.Nfft);
ovlap = 0;
SPECT_LTASS=[]; SPECT_LTASS_B=[];
for r = 1:size(speech,2)
   [SPECT_LTASS(:,r),freqs] = pwelch(speech(:,r),win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');SPECT_LTASS(:,r)=sqrt(SPECT_LTASS(:,r));
   [SPECT_LTASS_B(:,r),freqsB] = pwelch(speechB(:,r),win_,SYS.analysis_info.Nfft*ovlap,SYS.analysis_info.Nfft,SYS.signal_info.Fs,'power');SPECT_LTASS_B(:,r)=sqrt(SPECT_LTASS_B(:,r));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shape the noise to the speech being used
[~, SS, f_SS] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilt_SS( ...
    [], SYS.signal_info );

% Shape noise spectrum to match quiet zone leakage spectrum
[ ~, ~, QZS, QZS_low, f_QZS ] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_ZoneSpect( 'Quiet', ...
    [], SYS.Masker_Setup(1), SYS.Main_Setup, SYS.system_info, SYS.signal_info );

% Bright zone second order leakage spectrum
[ ~, ~, BZS, BZS_low, f_BZS ] ...
    = Broadband_Tools.Loudspeaker_Signal_Calculation.PreFilter_ZoneSpect( 'Quiet', ...
    [], SYS.Masker_Setup(1), SYS.Masker_Setup(1), SYS.system_info, SYS.signal_info );

if f_QZS ~= f_BZS, error('Frequency vectors don''t match'); end
lambda = 0; % lambda in {0,...,1} (Real)
ZS_rat1 = ( QZS.^(1-lambda) ) ./ ( BZS.^lambda );
ZS_rat_low1 = ( QZS_low.^(1-lambda) ) ./ ( BZS_low.^lambda );

lambda = 1/2; % lambda in {0,...,1} (Real)
ZS_rat2 = ( QZS.^(1-lambda) ) ./ ( BZS.^lambda );
ZS_rat_low2 = ( QZS_low.^(1-lambda) ) ./ ( BZS_low.^lambda );

lambda = 1; % lambda in {0,...,1} (Real)
ZS_rat3 = ( QZS.^(1-lambda) ) ./ ( BZS.^lambda );
ZS_rat_low3 = ( QZS_low.^(1-lambda) ) ./ ( BZS_low.^lambda );


% Adjust the noise to account for the aliasing caused by a limited number of loudspeakers
cheby_order = 9; Rp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(freqs - f_SS), error('Frequency vectors don''t match.');end

f_l = SYS.analysis_info.f_low;
f_c = Broadband_Tools.getAliasingFrequency(SYS.Masker_Setup(1))/2/pi*SYS.signal_info.c;
f_h = SYS.analysis_info.f_high;

f_band = (f_l<=freqs & freqs<=f_c);
f_band2 = (f_l<=f_QZS & f_QZS<=f_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
octSpace = SYS.signal_info.OctaveBandSpace;
P=[];PB=[];
for r = 1:size(speech,2)
    [P(:,r),f_]   = Tools.octaveBandMean(SPECT_LTASS(f_band,r),freqs(f_band),octSpace);
    [PB(:,r),f_B]  = Tools.octaveBandMean(SPECT_LTASS_B(f_band,r),freqsB(f_band),octSpace);
end
P = P./mean(P(:));
PB = PB./mean(PB(:));
[sp,f1] = Tools.octaveBandMean(SS(f_band),freqs(f_band),octSpace);
[q1,f2] = Tools.octaveBandMean(QZS(f_band2),f_QZS(f_band2),octSpace);
 q2     = Tools.octaveBandMean(QZS_low(f_band2),f_QZS(f_band2),octSpace);
 b1     = Tools.octaveBandMean(BZS(f_band2),f_QZS(f_band2),octSpace);
 b2     = Tools.octaveBandMean(BZS_low(f_band2),f_QZS(f_band2),octSpace);
 r11     = Tools.octaveBandMean(ZS_rat1(f_band2),f_QZS(f_band2),octSpace);
 r12     = Tools.octaveBandMean(ZS_rat_low1(f_band2),f_QZS(f_band2),octSpace);
 r21     = Tools.octaveBandMean(ZS_rat2(f_band2),f_QZS(f_band2),octSpace);
 r22     = Tools.octaveBandMean(ZS_rat_low2(f_band2),f_QZS(f_band2),octSpace);
 r31     = Tools.octaveBandMean(ZS_rat3(f_band2),f_QZS(f_band2),octSpace);
 r32     = Tools.octaveBandMean(ZS_rat_low3(f_band2),f_QZS(f_band2),octSpace);
 lp     = 1 ./ sqrt( 1 + (10^(Rp/10)-1)*chebyshevT(cheby_order,(f1/f_c)).^2);
if any(f1(1:end-1) - f2(1:end-1)), error('Frequency vectors don''t match.');end
f1([1 ])=[];
sp([1 ])=[];
q1([1 ])=[];
q2([1 ])=[];
b1([1 ])=[];
b2([1 ])=[];
r11([1 ])=[];
r12([1 ])=[];
r21([1 ])=[];
r22([1 ])=[];
r31([1 ])=[];
r32([1 ])=[];
lp([1 ])=[];
P([1 ],:)=[];
PB([1 ],:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wh = ones(size(f1)) * mean(P(:));
p = ones(size(f1)) * mean(P(:)) ./ (f1/mean(f1));
spr12lp = sp.*r12.*lp;
spr22lp = sp.*r22.*lp;
spr32lp = sp.*r32.*lp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [sp;q2;lp;wh;p;spr12lp;spr22lp;spr32lp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HQ=H;E=[]; Eps=[];H_final = HQ;
Eps_min = ones(1,size(HQ,1))*1e4;
E_min   = ones(size(P,2),size(HQ,1))*1e4;
for h = 1:size(HQ,1)
    HQ(h,:) = Tools.Correlated_Normalisation(mean(P,2),HQ(h,:));
    
    for g = db2mag(-30:.1:30)
        [E(:,h),Eps(h)]=Tools.COSHdist(HQ(h,:)*g,P);
        if Eps(h) <= Eps_min(h)
            H_final(h,:) = HQ(h,:)*g;
            E_min(:,h) = E(:,h);
            Eps_min(h) = Eps(h);
        end
    end
end
HB=H.*repmat(b1,size(H,1),1);EB=[]; EpsB=[];H_finalB = HB;
Eps_minB = ones(1,size(HB,1))*1e4;
E_minB   = ones(size(PB,2),size(HB,1))*1e4;
for h = 1:size(HB,1)
    HB(h,:) = Tools.Correlated_Normalisation(mean(PB,2),HB(h,:));
    
    for g = db2mag(-30:.1:30)
        [EB(:,h),EpsB(h)]=Tools.COSHdist(HB(h,:)*g,PB);
        if EpsB(h) <= Eps_minB(h)
            H_finalB(h,:) = HB(h,:)*g;
            E_minB(:,h) = EB(:,h);
            Eps_minB(h) = EpsB(h);
        end
    end
end

%%
Eps_min_dB = mag2db(Eps_min);
Eps_minB_dB = mag2db(Eps_minB);
Eps_mean_dB = mag2db(mean([Eps_min;Eps_minB]));
Eps_min_CI  = mag2db( exp( Tools.confidence_intervals(E_min,95,true) ) ).';
Eps_minB_CI = mag2db( exp( Tools.confidence_intervals(E_minB,95,true) ) ).';
Eps_mean_CI = mag2db( exp( Tools.confidence_intervals([E_min;E_minB],95,true) ) ).';

%%
if exist('fH'), if isvalid(fH), close(fH); end; end

colours = {[ ...            R G B  values
    0.2 0.2 1.0       ; ...       Bright Intelligibility
    1.0 0.0 0.0       ; ...       Quiet  Intelligibility
    0.6 0.0 0.6       ;];[ ...    Speech Intelligibility Contrast
    0.0 0.6 0.0       ];}; %      Bright Quality

LAB = rgb2lab(colours{1}(1:2,:));
[A,C] = cart2pol(LAB(:,2),LAB(:,3));
C = [74;84];
[a,b]=pol2cart(A,C);
LAB2 = [[60;40] a b];
colours{1}(1:2,:) = lab2rgb(LAB2);

% Figure Output Settings
DocumentPath = SYS.publication_info.DocumentPath;
print_fmt = 'pdf'; %figure image file format
print_res = 600; %dpi
plot_width = 88.9/10;% + 6.35/10 + 88.9/10; %IEEE full text width
aspect_ratio = 2/5;
FontSize = 8;
Font = 'Times';
lineWid = 0.5;
LegendHeightScaleFactor = 1.1;

% Latex Fonts
FF = SYS.publication_info.LaTeX_FontFamily; % Set fonts for latex interpreter
FFnums = SYS.publication_info.LaTeX_NumbersFontFamily; % Set number fonts for latex interpreter
FS = SYS.publication_info.FontSize;
latexFontSettings = ['{\fontfamily{' FF '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];
latexNumFontSettings = ['{\fontfamily{' FFnums '}' ...
    '\fontsize{' num2str(FS) 'pt}{' num2str(FS) 'pt}' ...
    '\selectfont '];

setLatexFont    = @(s) [latexFontSettings    s '}'];
setLatexNumFont = @(s) [latexNumFontSettings s '}'];


fH = figure(1);
fH.Color = 'w';

ax=gca; hold on;
yL = [-40 10];
xlim([0.5 5.5]);
ylim( yL );
ax.XTick = 1:5;
ax.YTick = -40:10:10;
xg = ax.XTick([2:end;2:end]) - 0.5 ; xg = xg(:).';
xgy = repmat([yL, flip(yL)]*1e3,1,(numel(ax.XTick)-1)/2);
plot( xg,xgy, '-','color',[0 0 0 0.5] );
yg = ax.YTick([1:end;1:end]); yg = yg(:).';
ygx = repmat([1 -1 -1 1]*1e3,1,ceil(numel(ax.YTick)/2)); ygx(numel(yg)+1:end)=[];
plot( ygx,yg, ':','color',[0 0 0 0.2] );


MarkSep = 0.2;
x = 1:5; x = x - MarkSep;
y = Eps_min_dB(end-4:end);
CI = Eps_min_CI(:,end-4:end);
errorbar(x,y,CI(1,:),CI(2,:),'.','color',colours{1}(1,:)); 

x = 1:5;
y = Eps_minB_dB(end-4:end);
CI = Eps_minB_CI(:,end-4:end);
errorbar(x,y,CI(1,:),CI(2,:),'.','color',colours{1}(2,:)); 

x = 1:5; x = x + MarkSep;
y = Eps_mean_dB(end-4:end);
CI = Eps_mean_CI(:,end-4:end);
errorbar(x,y,CI(1,:),CI(2,:),'.k');

ax.Box = 'on';
ax.TickDir = 'both';
ax.YMinorTick = 'on';
% ax.YGrid = 'on';
ax.TickLabelInterpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = setLatexFont('COSH Distance ($\mathrm{dB}$)');
ax.XTickLabel = {...
    '$\{\mathrm{ wh},\mathrm{lp}\}$'; ...
    '$\{\mathrm{  p},\mathrm{lp}\}$'; ...
    '\begin{tabular}{c} $\{\mathcal{IB},\mathrm{lp}\}$ \\ $\lambda{\grave{}}=0.0$\end{tabular}'; ...
    '\begin{tabular}{c} $\{\mathcal{IB},\mathrm{lp}\}$ \\ $\lambda{\grave{}}=0.5$\end{tabular}'; ...
    '\begin{tabular}{c} $\{\mathcal{IB},\mathrm{lp}\}$ \\ $\lambda{\grave{}}=1.0$\end{tabular}'; };
ax.YTickLabel = cellfun(@num2str,num2cell(ax.YTick(:)),'un',0);

ax.XTickLabel = cellfun(setLatexFont, ax.XTickLabel,'un',0);
ax.YTickLabel = cellfun(setLatexNumFont, ax.YTickLabel,'un',0);

grid off;
hold off;

ax.Units = 'centimeters';
ax.Position(3:4) = plot_width * [1 aspect_ratio];


ax2 = copyobj(ax,fH);
ax2.Color(4) = 0;
ax2.YLabel = [];
ax2.XTickLabel = [];
ax2.YTickLabel = [];
ax2.Children.delete;
ax.TickLength = [0 0];
ax2.XTick = ax.XTick(2:end) - 0.5;

drawnow; %pause(1.0);
tightfigadv;

fname = [SYS.publication_info.DocumentPath filesep 'COSH_Distances'];
print(fname,['-d' SYS.publication_info.print_fmt]);

%%
fprintf(['\nCOSH Distances in Decibels (dB)\n'...
    '\t\t\t'...
    '{sp}\t\t'...
    '{q''}\t\t'...
    '{lp}\t\t'...
    '{wh}\t\t'...
    '{p}\t\t\t'...
    '{IB(0),lp}\t'...
    '{IB(0.5),lp}\t'...
    '{IB(1),lp}\t'...
    '\n'...
    'Eps_COSH' '\t'     num2str(round(Eps_min_dB, 3,'significant')) '\n' ...
    'Eps_COSH_B' '\t'   num2str(round(Eps_minB_dB,3,'significant')) '\n' ...
    'Mean' '\t\t'       num2str(round(Eps_mean_dB,3,'significant')) '\n\n']);

%    
%%
Results = [Eps_min_dB, ...
           Eps_minB_dB, ...
           Eps_mean_dB];
newcom = '\\newcommand';
mac_pre = '\\';
mac_post = 'res';
mac_names = {'LTASS'; ...
             'QuietOne'; ...
             'LowPass'; ...
             'White'; ...
             'Pink'; ...
             'IBLambdaZero'; ...
             'IBLambdaHalf'; ...
             'IBLambdaOne'; };
mac_rownames = {''; ...
                'B'; ...
                'AVG';};
            
cols = numel(mac_names);
rows = numel(mac_rownames);
% mac_post = '';
% mac_names = num2cell(char((1:cols).' +64));
% mac_rownames = num2cell(char((1:rows).' +64));           
M = cellfun(@(a,b) [a,b],...
    repmat({ mac_pre   },cols,1),...
             mac_names ,...
    'UniformOutput',false);
macros = cellfun(@(a,b,c) [a,b,c],...
    repmat( M ,rows,1),...
    reshape(repmat( mac_rownames.' ,cols,1),[],1),...
    repmat({ mac_post  },rows*cols,1),...
    'UniformOutput',false);
      
FC1 = cellfun(@(m) [newcom '{' m '} {{' ],macros,'UniformOutput',false);
FC2 = [num2str(round(Results.',3,'significant')) repmat('}}\n',numel(macros),1)];

filecontent = cellfun(@(a,b) [a,b],FC1,mat2cell(FC2,ones(size(FC2,1),1),size(FC2,2)),'UniformOutput',false);

fid = fopen([SYS.publication_info.DocumentPath filesep SYS.publication_info.LatexMacrosFile],'w');
fprintf(fid,[filecontent{:}]);

fclose(fid);

Tools.MiKTeX_FNDB_Refresh;