
x = 0.01:0.01:20;

y = 3.64.*x.^(-0.8) - 6.5*exp(-0.6*(x-3.3).^2) + 10^-3*x.^4;

plot (x,y,'LineWidth',3); set(gca,'XScale','log'); grid on;

hold on;
x = [1.25 1.25];
plot (x,[-20 180],'black','LineWidth',3);

x = [8 8];
plot (x,[-20 180],'black', 'LineWidth',3);

hold off;