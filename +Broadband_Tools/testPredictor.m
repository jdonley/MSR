
i=50;
figure(1)
plot(B(:,i+4));
hold on;
plot(y(:,i));
hold off;

figure(2)
i2=i+1;
plot(B(:,i2+4));
hold on;
plot(y(:,i2));
hold off;

%%
h = a(i,:);
[z,p,k] = tf2zpk(h,1);
hmin = polystab(h)*norm(h)/norm(polystab(h));
a2 = zp2tf(z,p,k);
plot(a(i,:));hold on;
plot(hmin);hold on;
 plot(a2);hold off;
%%
a(isnan(a))=0;
z=zeros(size(a,1),size(a,2)-1);
for k = 1:size(a,1)
z_tmp = tf2zpk(a(k,:),1);    
if ~isempty(z_tmp)
z(k,:) = z_tmp;
end
end


%%
[z,p,k] = tf2zpk(a(i,:),1);
figure(3)
plot(abs(z))

z2 = z;
z2_unstab = z2(abs(z2)>1);
z2(abs(z2)>1) = (2-abs(z2_unstab)) .* exp(1i*angle(z2_unstab));
hold on;
plot(abs(z2));

figure(4)
scatter(real(z),imag(z)); hold on;
scatter(real(z2_unstab),imag(z2_unstab)); hold on
scatter(real(z2),imag(z2)); hold off


%%

Nfft=size(y,1);
window_ = [hanning(Nfft-1);0];
Y = overlapadd(F_p((end-m+1):end,:).',window_,Nfft/2);
Yorig = overlapadd(B(1:end/2,5:end).',window_,Nfft/2);

Yorig(end+1:numel(Y))=0;
Y(isnan(Y))=0;

p = rms(Yorig);
Yorig=Yorig/p;
Y=Y/p;
%
s = 0:0.01:2;
Prms = [];
for s_ = 1:numel(s)
scale = s(s_);
Prms(s_) = rms(Yorig-Y*scale);
end

% scale = s(find(Prms==min(Prms)))
figure(3)
plot(Yorig,'k'); hold on;
plot(Y*scale,'r'); hold on

plot(Yorig-Y*scale,'b'); hold on
hold off;

mag2db(Prms(find(Prms==min(Prms))))
