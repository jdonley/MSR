function [ y_ZS, y_ZS_low, spect, spect_pass, frqs ] = PreFilter_ZoneSpect( Zone, y, setup, multizone_setup, system_info, signal_info )
%PREFILTER_QZS Shape signal to the Quiet Zone Spectrum

%obtain the quiet zone spectrum from the database of reproduction estimates
[spect,frqs,DBfrqs] = Soundfield_Database.getZoneSpectrumFromDB( Zone, multizone_setup, system_info.LUT_resolution, system_info.Drive, signal_info.Fs );
spect_pass = spect;

f_cutoff = Broadband_Tools.getAliasingFrequency(setup) * signal_info.c / (2*pi);
%remove sharp increase from aliasing
if any(frqs>f_cutoff)
    spect_high=spect(frqs>f_cutoff);
    spect_pass(frqs>f_cutoff) = spect_high(1);
end

%limit to bandwidth of interest
f_low = max([signal_info.f_low,min(DBfrqs)]);
if any(frqs<f_low)
    spect_low=spect(frqs<f_low);
    spect_pass(frqs<f_low) = spect_low(end);
end

%filter the signal
y_ZS     = Tools.ArbitraryOctaveFilt(y,     spect, frqs, signal_info.Nfft, signal_info.Fs, signal_info.OctaveBandSpace);
y_ZS_low = Tools.ArbitraryOctaveFilt(y, spect_pass, frqs, signal_info.Nfft, signal_info.Fs, signal_info.OctaveBandSpace);

end

