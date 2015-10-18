f = logspace(log10(20),log10(12500),1000);

thresh1 = Perceptual_Tools.Threshold_in_Quiet(f,'approx');
thresh2 = Perceptual_Tools.Threshold_in_Quiet(f,'ISO226');
thresh3 = Perceptual_Tools.Loudness(f,0);

plot(f,thresh1);hold on;
plot(f,thresh2);
plot(f,thresh3);hold off;

set(gca,'XScale','log'); legend show;
grid