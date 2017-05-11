diameter_	= [0.5:7]*2.54/100;		% convert inches to meters
Xmax_		= [0.5:7]/1000;		% convert to meters
[diameter,Xmax]=meshgrid(diameter_,Xmax_);

drivers = 1; % number of drivers
C  = 345;			% speed of sound m/s
Ro = 1.18;
SPL = 100; % dB SPL @ 1m

Hz = ((10^((SPL - 112) / 10)) ./ (pi^5 * Ro / 4 / C * (drivers * Xmax).^2 .* (0.83 * diameter).^4)).^0.25

round(Hz)