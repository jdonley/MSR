clc;clear;

%%
filepath = ['Z:\+Recordings\' ...
    '+Physical_World\' ...
    '+CIRCLEarray_1.3mPerpDist_180degCentre\' ...
    '+24Genelec8010ASpkrs_4.2616mLen\' ...
    '+0Bx_0.6By_0Qx_-0.6Qy\' ...
    '+VSrc_pw_0deg_0m\' ...
    '+Database_512f_32w\' ...
    '+10x10x10Dim_1Ab\' ...
    '32Rec_5x5x5Ctr\' ...
    '+150Hz-8000Hz_-InfdB_10000weight__method_ZoneWeightMaskerAliasCtrl\'];

% filepath = ['Z:\+Recordings\' ...
%     '+CIRCLEarray_1.3mPerpDist_180degCentre\' ...
%     '+24Genelec8010ASpkrs_4.2616mLen\' ...
%     '+0Bx_0.6By_0Qx_-0.6Qy\' ...
%     '+VSrc_pw_0deg_0m\' ...
%     '+Database_512f_32w\' ...
%     '+10x10x10Dim_1Ab\' ...
%     '32Rec_5x5x5Ctr\' ...
%     '+150Hz-8000Hz_-InfdB_10000weight__method_NoMask\'];
%%
extractWithText = 'Original';
files = Tools.getAllFiles(filepath);
filenames={};
for f = 1:length(files)
    if ~isempty(findstr(files{f},'Sweep'))
        continue;
    end
    if ~isempty(findstr(files{f},extractWithText))
        [~,tmp] = fileparts(files{f});
        tmp = flip(tmp);
        tmp = flip(tmp(length(extractWithText)+2:end));
        filenames{end+1,1}=tmp;
    end
end

%%
pesqVal = [];
for f = 1:length(filenames)
    audio_filename = filenames{f};
    
    [orig,fs] = audioread([filepath audio_filename '_Original.WAV']);
    load([filepath audio_filename '_Bright.mat'])
    ref = Rec_Sigs_B;
    %[ref,fs] = audioread([filepath audio_filename '_Reference.WAV']);
    
    if size(ref,1)<size(ref,2), ref=permute(ref(1,:),[2 1]); end
    
    ref_norm = Broadband_Tools.power_norm(double(orig),double(ref),fs,[150 7500]);
    orig_deci=decimate(orig,fs/16000);
    ref_norm_deci=decimate(ref_norm,fs/16000);
    [b,a]=butter(6,[150 7500]./(16000/2));
    orig_deci_filt = filter(b,a,orig_deci);
    ref_norm_deci_filt = filter(b,a,ref_norm_deci);
    pesqVal(f) = Tools.pesq_mex_vec(orig_deci_filt,ref_norm_deci_filt,16000,1);

end

meanPESQMOSval = mean(pesqVal);

fprintf('\nMean PESQ MOS value: %.2f\n\n',meanPESQMOSval);