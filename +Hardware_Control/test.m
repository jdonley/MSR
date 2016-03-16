clc;
clear;
playrec('reset');

%% Settings
dev_model = 'ASIO Hammerfall DSP';
fs = 48000;
Nspkrs = 24;

%% Find and initialise hardware
devs = playrec('getDevices');
for d = 1:length(devs)
    if strcmpi(devs(d).name,dev_model)
        dev = devs(d);
        break;
    end
end

playrec( 'init', fs, dev.deviceID, dev.deviceID );

%% Use
test_duration = 0.5;%seconds
y = randn( test_duration*fs, 1 );
y = y ./ max(abs(y(:)));
y_multi = reshape( [repmat( [y; zeros(test_duration*fs * Nspkrs,1)], Nspkrs - 1, 1); y], ...
            test_duration * fs * Nspkrs, ...
            Nspkrs);
chans = 1:Nspkrs;
playrec('play',y_multi,chans);
