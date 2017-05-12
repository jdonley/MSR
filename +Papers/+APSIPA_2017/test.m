diameter	= linspace(1,10,800)/100;		% convert inches to meters
% Xmax_		= linspace(0.0,2,800)/100;		% convert to meters
Xmax = diameter * 10/100; 
% [diameter,Xmax]=meshgrid(diameter_,Xmax_);

drivers = 1; % number of drivers
C  = 345;			% speed of sound m/s
Ro = 1.18;
SPL = 100; % dB SPL @ 1m

Hz = ((10^((SPL - 112) / 10)) ./ (pi^5 * Ro / 4 / C * (drivers * Xmax).^2 .* (0.83 * diameter).^4)).^0.25;

round(Hz);

% 
% surf(diameter_*100, Xmax_*100, Hz/1e3,'linestyle','none')
% set(gca,'ZScale','log');
% xlabel('Diameter')
% ylabel('Max Excursion')
% zlabel('Frequency (kHz)')

plot( Hz/1e3, diameter*100 )
ylabel('Diameter')
xlabel('Frequency (kHz)')
 set(gca,'XScale','log');
 xlim([0.05 8]);

