function IRgains = getGainAdjustments( IRs, f_band, fs )

% Gain compensation
Nspkrs = size(IRs,2);
pkIRval = max(abs(IRs));
[~,loudestSpkr] = max(pkIRval);

IRgains = zeros(1,Nspkrs);
for spkr = 1:Nspkrs
    [~,IRgains(spkr)] = Broadband_Tools.power_norm( ...
        IRs(:,loudestSpkr), ...
        IRs(:,spkr), ...
        fs, ...
        f_band);
end

end