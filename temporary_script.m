tic;
SYS = Current_Systems.IEEELetters2017_System_B;

subSYS = SYS;
subSYS.Main_Setup(1)=[];
subSYS.Room_Setup(1)=[];
subSYS.signal_info.method = 'Clean';
MicSigPath = Broadband_Tools.getMicrophoneSignalPath(subSYS);
MicSigFiles = Tools.keepFilesFromFolder( Tools.getAllFiles(MicSigPath), subSYS.signal_info.speech_filepath);
MicSigFiles(~contains(MicSigFiles,'.mat'))=[];
data = load(MicSigFiles{:});
fs = data.fs;

SNRdb = 60;

room = SYS.Room_Setup(2);
Q = room.NoReceivers;
c = SYS.signal_info.c;
fmax = 2000;
R = SYS.Main_Setup(1).Radius;
Nmax = ceil(2*pi*fmax/c*R);

WIN = sqrt(hann(512,'p'));

S=[];mic_sigs=[];
for l = 1:size(data.mic_signals,2)

    mic_sigs(:,l) = data.mic_signals(:,l);%awgn(data.mic_signals(:,l).',SNRdb,'measured').';
    [S(:,:,l),ff,tt] = spectrogram(mic_sigs(:,l),WIN,256,512,fs);
    ss(:,:,l) = enframe(mic_sigs(:,l),WIN,256);
    SStmp = fft(permute(ss(:,:,l),[2 1 3]));
    SS(:,:,l) = SStmp(1:end/2+1,:);
    
end
%% Determine soundfield coefficients for mics that aren't colocated 

% mLocs = room.ReceiverPositions - repmat(room.Reproduction_Centre([2 1 3]),Q,1);
% [th,r] = cart2pol( mLocs(:,1), mLocs(:,2) );
% 
% B=zeros(numel(ff),numel(tt),2*Nmax+1);
% badFreq=[];
% fprintf('\t Completion: '); Tools.showTimeToCompletion; startTime=tic;
% for t_ = 301%1:numel(tt)
%     t = tt(t_);
%     
%     for f_ = 17%2:numel(ff)
%         f = ff(f_);        
%         k = 2*pi*f/c;
%         N = ceil(k*R);
%         
%         Y = squeeze(S(f_,t_,:)).';
%         Y = repmat(Y,2*N+1,1);
%         [rr,NN] = meshgrid( r,-N:N);
%         thth    = meshgrid(th+pi/2,-N:N);
%         
%         H = 1i/4*besselh(NN, k * rr) ;
%         T = 1/Q ./ H .* exp( -1i*NN.*thth );
%         
%         
%         BB = sum(Y .* T, 2);
%         B(f_,t_,:) = [zeros(Nmax-N,1);BB;zeros(Nmax-N,1)];
%         
%         
%         Tools.showTimeToCompletion(f_/numel(ff), [], [], startTime );
%     end
% end
% 
% 
% toc;
% % Determine Loudspeaker Weights
% setup = SYS.Main_Setup(1);
% SpkrLocs = setup.Loudspeaker_Locations;
% L = size(SpkrLocs,1);
% phi = setup.Speaker_Arc_Angle / 180 * pi;
% L_ = L/2;
% delta_phi_s = phi / L_;
% 
% [SpkrLocs(:,1),SpkrLocs(:,2)]=pol2cart(SpkrLocs(:,1),SpkrLocs(:,2));
% 
% [th,r] = cart2pol( SpkrLocs(:,1), SpkrLocs(:,2) );
% 
% 
% ff(ff>2500)=[];
% c = SYS.signal_info.c;
% 
% 
% 
% X=[]; 
% for t_ = 301%:numel(tt)
%     for f_ = 17%2:numel(ff)
%         f = ff(f_);
%         
%                
%         k = 2*pi*f/c;
%         N = ceil(k*R);
%         
%         [rr,NN] = meshgrid( r,-N:N);
%         thth    = meshgrid(th,-N:N);
%         
%         beta = repmat( squeeze(B(f_,t_,:)), 1, L );
% 
%         
%         X(f_,:) = sum(  beta((-N:N)+Nmax+1,:) .* exp(1i*NN.*thth) * delta_phi_s ,1 );
%         
%         
%     end
%     disp(t_)
% end


%% Determine position or replicated soundfield 
searchFieldRes = 30;
x = 0.2 : 1/searchFieldRes : room.Room_Size(2);
y = 0.2 : 1/searchFieldRes : room.Room_Size(1);

setup = SYS.Main_Setup(1);
c = SYS.signal_info.c;
SpkrLocs = setup.Loudspeaker_Locations;
L = size(SpkrLocs,1);
[SpkrLocs(:,1),SpkrLocs(:,2)]=pol2cart(SpkrLocs(:,1),SpkrLocs(:,2));
 SpkrLocs = [SpkrLocs zeros(L,1)] + repmat(room.Reproduction_Centre([2 1 3]),size(SpkrLocs,1),1);
                
 [xx_,yy_] = meshgrid( x , y );
 
 %%%
 
 SLx = repmat(permute(SpkrLocs(:,1),[2 3 1]),[size(xx_) 1]);
 SLy = repmat(permute(SpkrLocs(:,2),[2 3 1]),[size(yy_) 1]);
 xx = repmat(xx_,1,1,L) - SLx;
 yy = repmat(yy_,1,1,L) - SLy;
 
 [th,r] = cart2pol(xx,yy);
 
%  FF = (ff(end)*((2:numel(ff)).'-1)/numel(ff));
 FF = ff(2:end);

 ffI = FF<1000 | FF>2500;
FF(ffI)=[];
 k = 2*pi*FF/c;
 kk = repmat( permute(k,[2 3 4 1]),[size(r) 1]);
 rr = repmat(r,[1 1 1 numel(FF)]);
 
%  H = 1i/4*besselh(0,kk.*rr); %2D
 H = exp(1i*kk.*rr) ./ (4*pi*rr); %3D;
 exps = exp( 1i*c*kk ); 
 
 %%%
 %%
 VSlocpol = SYS.Main_Setup(2).Loudspeaker_Locations;
 [VSloc(1), VSloc(2)] = pol2cart(VSlocpol(1),VSlocpol(2));
 VS = SYS.Room_Setup(2).Reproduction_Centre([2 1]) + VSloc;
 
 [~,c_true]=min(abs(x-VS(1)));
 [~,r_true]=min(abs(y-VS(2)));

 
 k2 = 2*pi*ff/c;
 kk2 = repmat(permute(k2,[2 1]),[size(r,3) 1]);
 rr2 = repmat(squeeze(r(r_true,c_true,:)),[1 numel(ff)]);
 
%  H2 = 1i/4*besselh(0,kk2.*rr2); %2D
 H2 = exp(1i*kk2.*rr2) ./ (4*pi*rr2); %3D
 exps2 = exp( 1i*c*kk2 );
 
 
 
 ignoreLevel = max(abs(S(:)))*db2mag(-20);
 Ssrc = [];
 FIELD=[];FIELDTOT=[]; tic;
 for t_ = 1:1:numel(tt)
     
%      if max(max(abs(S(~ffI,t_,:)))) < ignoreLevel % Ignore speech level that is lower than -20dB of peak value
%          continue
%      end
%      XX = repmat( permute(S(~ffI,t_,:),[4 2 3 1]), [size(xx_) 1 1 ] );
     XX2 = squeeze(S(:,t_,:)).';
     
     
     VirtualSenseSig = XX2 ./ conj(H2);% .* exps2;
     
     VirtualSenseSig(:,1) = 0;
     VirtualSenseSig(:,end) = real(VirtualSenseSig(:,end));
     
     Ssrc(:,t_) = sum(VirtualSenseSig);
     
% %      FLD = sum( XX .* H .* exps, 3 );
% %      FIELD(:,:,t_) = (sum( abs(FLD), 4 )).^2; % Squaring the absolute sum helps when noise is present
%      FLD = sum( XX ./ conj(H), 3 );
%      FIELD(:,:,t_) = (sum( abs(FLD), 4 )); % Squaring the absolute sum helps when noise is present
%      
%       disp(t_)
%      
%      avgingLen = ceil( SYS.signal_info.rir_duration /diff(tt(1:2)))+1;
%      if t_-avgingLen > 0
%          FIELDTOT = (squeeze(sum(FIELD(:,:,t_-avgingLen:t_),3))); % Using a longer time segment helps when reverberation exists
%      else
%          FIELDTOT = (squeeze(sum(FIELD(:,:,1:t_),3)));
%      end
% %      FIELDTOT = (FIELD(:,:,t_));
%      
%      
%      % Search for point source origin via soundfield correlations
%      figure(2); hold off;
%      
%      [fldMax,I]=max(FIELDTOT(:));
%      [r_,c_]=ind2sub(size(FIELDTOT),I);
%      scatter3(VS(1),VS(2),fldMax,'ro'); hold on;
%      scatter3(x(c_),y(r_),fldMax,'k.'); hold on;
%      
%      surf(x,y,FIELDTOT,'lines','no');
%      view(2);axis equal;
%      xlim([min(x) max(x)]);ylim([min(y) max(y)]);
%      drawnow;
 end
 
 
 toc;
 
 Ssrcfft = [Ssrc; conj(Ssrc(end-1:-1:2,:))];
 
 srcSigFrm = ifft(Ssrcfft);
 
  srcSig = Broadband_Tools.OverlapAdd((srcSigFrm .* Tools.repmatmatch(WIN,srcSigFrm)).',0.5);
%  srcSig = overlapadd(srcSigFrm.',WIN,256);
 
%  [b,a] = cheby1(6,1,[250 2500]/8000);
%  srcSig = filter(b,a,srcSig);
 
 load('M:\MSR\+Miscellaneous\+TestAudio_Files_InvFilts\inverseESS.mat');
 imp = Tools.extractIR(srcSig,invY);
 plot(imp)