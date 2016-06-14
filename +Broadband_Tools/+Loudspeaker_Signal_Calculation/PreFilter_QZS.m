function [ y_QZS, y_QZS_low, spect, spect_low, frqs ] = PreFilter_QZS( y, setup, multizone_setup, system_info, signal_info )
%PREFILTER_QZS Shape signal to the Quiet Zone Spectrum

[spect,frqs] = Soundfield_Database.getQuietZoneSpectrumFromDB( multizone_setup, system_info.LUT_resolution, system_info.Drive, signal_info.Fs );
spect_low = spect;

f_cutoff = Broadband_Tools.getAliasingFrequency(setup) * signal_info.c / (2*pi);
%remove sharp increase from aliasing
spect_high=spect(frqs>f_cutoff);
spect_low(frqs>f_cutoff) = spect_high(1);

y_QZS     = Tools.ArbitraryOctaveFilt(y,     spect, frqs, signal_info.Nfft, signal_info.Fs, 1/6);
y_QZS_low = Tools.ArbitraryOctaveFilt(y, spect_low, frqs, signal_info.Nfft, signal_info.Fs, 1/6);

end

