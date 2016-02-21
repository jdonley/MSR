function saveLoudspeakerSignals( path, name, ext, spkr_signals, spkr_count, original, orig_name, orig_ext, fs )
%SAVELOUDSPEAKERSIGNALS Summary of this function goes here
%   Detailed explanation goes here
if isempty(spkr_count)
    spkr_count = size(spkr_signals,2);
end

filenumbers = num2str((1:spkr_count)');
filenumbers(filenumbers==' ') = '_';
fullpath = [repmat([path name], [spkr_count 1]) ...
    filenumbers ...
    repmat(ext, [spkr_count 1]) ];

if ~exist(path,'dir'); mkdir(path); end

for spkr = 1:spkr_count
    audiowrite(fullpath(spkr,:), spkr_signals(:, spkr), fs);
end

if ~isempty(original)
    audiowrite([path ...
        orig_name '_' 'Original' ...
        orig_ext], original, fs);
end

end

