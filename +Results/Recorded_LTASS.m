%% Load current system
SYS = Current_Systems.loadCurrentSRsystem;

%%
recPath = Results.getRecordingsPath( SYS );
% Hardware_Control.getRealRecordingsPath
%% Get all audio
files = Tools.getAllFiles( recPath );
speechB=[];
speechQ=[];
fs=[];
F = length(files);
for file = 1:F
    if strfind(lower(files{file}),'bright.mat')
        spB = load(files{file});
        if isfield(spB,'fs'), fs=spB.fs; else
        fs = SYS.signal_info.Fs; end;
        spB_RMS = mean(rms(spB.Rec_Sigs_B,2));
        sp = spB.Rec_Sigs_B ./ spB_RMS;
        speechB = [speechB, sp];
    elseif strfind(lower(files{file}),'quiet.mat')
        spQ = load(files{file});
        sp = spQ.Rec_Sigs_Q ./ spB_RMS;
        speechQ = [speechQ, sp];
    end
end

%%
[spectB,frqsB] = Tools.LTASS(speechB.',SYS.signal_info.Nfft,fs);
[spectQ,frqsQ] = Tools.LTASS(speechQ.',SYS.signal_info.Nfft,fs);


%%
frqAlias = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup)/2/pi*SYS.signal_info.c;
frqAliasOld = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup,'old')/2/pi*SYS.signal_info.c;

%%
spQ_RMS = mean(rms(spectQ(frqsQ < frqAlias,:),2));

%%
for m = 1:2
S = SYS;
S.Main_Setup = SYS.Masker_Setup;
S.signal_info.L_noise_mask = 0;
S.signal_info.method = S.signal_info.methods_list{S.signal_info.methods_list_masker(m)};

recPath = Results.getRecordingsPath( S );
% Hardware_Control.getRealRecordingsPath

files = flip(Tools.getAllFiles( recPath ));
maskerB=[];
maskerQ=[];
fs=[];
F = length(files);
for file = 1:F
    if strfind(lower(files{file}),'bright.mat')
        mskB = load(files{file});
        if isfield(mskB,'fs'), fs=mskB.fs; else
        fs = SYS.signal_info.Fs; end;
        msk = mskB.Rec_Sigs_B ./ mskQ_RMS;
        maskerB = [maskerB, msk];
    elseif strfind(lower(files{file}),'quiet.mat')
        mskQ = load(files{file});
        mskQ_RMS = mean(rms(mskQ.Rec_Sigs_Q,2));
        msk = mskQ.Rec_Sigs_Q ./ mskQ_RMS;
        maskerQ = [maskerQ, msk];
    end
end
[spectMskB(:,:,m),frqsMskB(:,:,m)] = Tools.LTASS(maskerB.',SYS.signal_info.Nfft,fs);
[spectMskQ(:,:,m),frqsMskQ(:,:,m)] = Tools.LTASS(maskerQ.',SYS.signal_info.Nfft,fs);
spectMskB(:,:,m) = spectMskB(:,:,m) ./ mean(rms(spectMskQ(frqsMskQ(:,1,m) < frqAlias,:,m),2)) .* spQ_RMS;
spectMskQ(:,:,m) = spectMskQ(:,:,m) ./ mean(rms(spectMskQ(frqsMskQ(:,1,m) < frqAlias,:,m),2)) .* spQ_RMS;
end


%%
figure(1);
plB = plot(repmat(frqsB,1,size(spectB,2)),mag2db((spectB))); hold on;
plB = plot([1 1]*frqAlias,[-160 -20]); hold on;
plB = plot([1 1]*frqAliasOld,[-160 -20]); hold on;
plQ = plot(repmat(frqsQ,1,size(spectQ,2)),mag2db(mean(spectQ,2))); hold on;
plQ = plot(repmat(frqsMskQ(:,:,1),1,size(spectMskQ(:,:,1),2)),mag2db(mean(spectMskQ(:,:,1),2)),'LineWidth',2); hold on;
plQ = plot(repmat(frqsMskQ(:,:,2),1,size(spectMskQ(:,:,2),2)),mag2db(mean(spectMskQ(:,:,2),2)),'LineWidth',2); hold on;
hold off;
set(gca,'XScale','log');
xlim([0.1 10]*1000);
ylim([-160 -20]);
grid on; grid minor;

%% Find some distances 

D_COSH = zeros(size(spectMskQ,2),2);
for m = 1:2
    for r = 1:size(spectMskQ,2)
       D_COSH(r,m) = distchpf(spectQ(frqsQ < frqAlias,r).',spectMskQ(frqsMskQ(:,1,m) < frqAlias,r,m).'); 
    end
end
D_COSH
mean(D_COSH,1)