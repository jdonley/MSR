function [ sig_matched, adjustVal ] = power_norm( sig_orig, sig_diff, Fs, band)
%POWER_NORM Summary of this function goes here
%   Detailed explanation goes here
if length(band) == 1
    % Calculate octave band
    f_centre = band;
    fd = (2^0.5);
    band = f_centre .* [1/fd, fd];
end

bpow_in   = mean(bandpower( sig_orig, Fs, band));
bpow_diff = mean(bandpower( sig_diff, Fs, band));
adjustVal = mean( db2mag( pow2db( bpow_in ./ bpow_diff ) ) );
sig_matched = sig_diff * adjustVal;

end

