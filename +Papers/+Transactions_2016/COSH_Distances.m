clear;
%%
SYS = Current_Systems.loadCurrentSRsystem;
res_path = Results.getRecordingsPath(SYS);

files = Tools.getAllFiles( res_path );
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
if any(f1 - f2), error('Frequency vectors don''t match.');end
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
    'Mean' '\t\t'       num2str(round(mean([Eps_min_dB;Eps_minB_dB]),3,'significant')) '\n\n']);

%    
%%
Results = [Eps_min_dB, ...
           Eps_minB_dB, ...
           mean([Eps_min_dB;Eps_minB_dB])];
newcom = '\\newcommand';
macros = {'\\LTASSres'; ...
          '\\QuietOneres'; ...
          '\\LowPassres'; ...
          '\\Whiteres'; ...
          '\\Pinkres'; ...
          '\\LTASSQuietOneLowPassres'; ...
          '\\IBLambdaHalfres'; ...
          '\\IBLambdaOneres';...
          '\\LTASSBres'; ...
          '\\QuietOneBres'; ...
          '\\LowPassBres'; ...
          '\\WhiteBres'; ...
          '\\PinkBres'; ...
          '\\LTASSQuietOneLowPassBres'; ...
          '\\IBLambdaHalfBres'; ...
          '\\IBLambdaOneBres';...
          '\\LTASSAVGres'; ...
          '\\QuietOneAVGres'; ...
          '\\LowPassAVGres'; ...
          '\\WhiteAVGres'; ...
          '\\PinkAVGres'; ...
          '\\LTASSQuietOneLowPassAVGres'; ...
          '\\IBLambdaHalfAVGres'; ...
          '\\IBLambdaOneAVGres';};
      
FC1 = cellfun(@(m) [newcom '{' m '} {{' ],macros,'UniformOutput',false);
FC2 = [num2str(round(Results.',3,'significant')) repmat('}}\n',numel(macros),1)];

filecontent = cellfun(@(a,b) [a,b],FC1,mat2cell(FC2,ones(size(FC2,1),1),size(FC2,2)),'UniformOutput',false);

fid = fopen([SYS.publication_info.DocumentPath filesep SYS.publication_info.LatexMacrosFile],'w');
fprintf(fid,[filecontent{:}]);

fclose(fid);

Tools.MiKTeX_FNDB_Refresh;