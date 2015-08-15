[Amplitudes, Phases, Frequencies] = Broadband_Tools.FFT_custom('+Miscellaneous\Test16k.wav', 1024, 16000);

surf(Frequencies, 1:size(Amplitudes,1), Amplitudes, 'LineStyle', 'none')
ylim([0 size(Amplitudes,1)])
set(gca, 'XScale', 'lin')
view(2)
caxis([0 0.01])