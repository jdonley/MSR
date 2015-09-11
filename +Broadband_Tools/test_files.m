% %%
% freqs=logspace(log10(150),log10(8000),16);
% t=0:1/16000:0.5-1/16000;
% y=zeros(16,length(t));
% for i=1:16
% y(i,:)=sin(2*pi*freqs(i)*t);
% end
% ze=zeros(size(y));
% res=[y ze]';
% res=res(:);
% audiowrite('sinusoids_150-8000Hz.WAV',res,16000);
% 
% 
% %%
% close all;
% 
% fs=16000;
% len=size(Rec_Sigs_B_500,2)-1;
% NFFT = len;
% f=fs*(1:NFFT/2-1)/NFFT;
% 
%  hold on;
%  for i=1:32
%      y=Rec_Sigs_B_500(i,:);     
%      Y=fftshift(fft(y,NFFT));
%      plot(f,abs(Y(end/2+1:end-1)));
% %      [pk,I]=max(abs(Y(end/2+1:end-1)));
% %      plot(f(I),pk,'+');
%      
%      y=Rec_Sigs_B_1000(i,:);
%      Y=fftshift(fft(y,NFFT));
%      plot(f,abs(Y(end/2+1:end-1)));
% %      [pk,I]=max(abs(Y(end/2+1:end-1)));
% %      plot(f(I),pk,'+');
%      
%      y=Rec_Sigs_B_4000(i,:);
%      Y=fftshift(fft(y,NFFT));
%      plot(f,abs(Y(end/2+1:end-1)));
% %      [pk,I]=max(abs(Y(end/2+1:end-1)));
% %      plot(f(I),pk,'+');
%      
%      set(gca,'XScale','log');
%      xlim([100 10000]);
%  end
%  hold off;

 %% chirp
 close all;
 figure(2)
 
 fs=16000;
len=length(Rec_Sigs_B(1,:))-1;
NFFT = len;
f=fs*(1:NFFT/2-1)/NFFT;

 hold on;
 for i=1:32

%      y=data;
%      Y=fftshift(fft(y,NFFT));
%      Z=abs(Y(end/2+1:end-1));
%      plot(f,Z);
%        [pk,I]=findpeaks(Z,f,'NPeaks',15,'MinPeakHeight',200,'MinPeakDistance',40);
%        plot(I,pk);


     y=Rec_Sigs_B(i,:);     
     Y=fftshift(fft(y,NFFT));
     Z=abs(Y(end/2+1:end-1));
     plot(f,Z);
      [pk,I]=findpeaks(Z,f,'NPeaks',15,'MinPeakHeight',500,'MinPeakDistance',40);
      plot(I,pk);

     set(gca,'XScale','log');
     box on;
     grid on;
     xlim([100 10000]);
 end
 hold off;
 