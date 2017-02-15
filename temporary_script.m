f=1000;
c=343;
k = 2*pi*f/c;

fmid=1000;
kset=2*pi*fmid/c;

d=0.02; %kd<1

nT=5;
x = linspace(-c/f*nT,c/f*nT,500);
x1 = x-d/2;
x2 = x+d/2;
[xx,yy]=meshgrid(x);
xx1=meshgrid(x1,x);
xx2=meshgrid(x2,x);
[tt,rr]=cart2pol(xx,yy);
[tt1,rr1]=cart2pol(xx1,yy);
[tt2,rr2]=cart2pol(xx2,yy);

% r = abs(x);
% r2 = abs(x2);


a1 = besselh(0,k*rr);

% compshift = exp(1i*pi*(1/2-2*d*f/c));
% b1 = besselh(0,k*rr) * compshift;
% c1 = besselh(0,k*rr2)*exp(1i*(1+2*d*f/c)*pi) * compshift;

b1 = besselh(0,k*rr1) * exp( -1i*pi/2) * exp( 1i*pi*2*d*f/c);
c1 = besselh(0,k*rr2) * exp( 1i*pi/2) ;

b1 = besselh(0,k*rr1)  * (-1i/(k*d*2)) * exp( 1i*d/2*k);
c1 = besselh(0,k*rr2)  * (1i/(k*d*2)) * exp( -1i*d/2*k);


c1 = besselh(0,k*rr2)  * exp( 1i*(k*d - pi)/2) / (2*k*d);
b1 = besselh(0,k*rr1)  * exp( 1i*(pi - k*d)/2) / (2*k*d);


figure(11)
surf( real(a1) ,'linestyle','none');view(2)
figure(1)
surf( real(b1+c1) ,'linestyle','none');view(2)

figure(2)
plot(x,real(a1(end/2,:)),'k'); hold on;
% plot(x,real(-a1(end/2,:)),':k'); hold on;
plot(x,real(b1(end/2,:)),'--r');hold on;
plot(x,real(c1(end/2,:)),'--m');hold on;
plot(x,real(b1(end/2,:)+c1(end/2,:)),'b');
grid on; grid minor;
hold off
axis([-c/f*nT c/f*nT -2 2])