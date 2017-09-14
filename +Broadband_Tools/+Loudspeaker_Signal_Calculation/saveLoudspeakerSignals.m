function saveLoudspeakerSignals( path, name, ext, spkr_signals, spkr_count, original, orig_name, orig_ext, fs )
%SAVELOUDSPEAKERSIGNALS Summary of this function goes here
%   Detailed explanation goes here
if isempty(spkr_count)
    spkr_count = size(spkr_signals,2);
end

% filenumbers = num2str((1:spkr_count)');
% filenumbers(filenumbers==' ') = '_';
% fullpath = [repmat([path name], [spkr_count 1]) ...
%     filenumbers ...
%     repmat(ext, [spkr_count 1]) ];

if ~exist(path,'dir'); mkdir(path); end

% for spkr = 1:spkr_count
%     audiowrite(fullpath(spkr,:), spkr_signals(:, spkr), fs, 'BitsPerSample',64);
% end

fullpath_multiCh = [path name num2str(spkr_count) 'Ch' ext];

% MATLABs audiowrite fails on high channel counts (>255) so we use
% VOICEBOX's writewav() and readwav() functions
% audiowrite(fullpath_multiCh, spkr_signals, fs, 'BitsPerSample',64);
writewavOpts = ['V' ...% 64 bit floating point
                'g' ...% save scaling value to the wave file to restore scale upon readwav()
                'z'];  % No PCM offset correction
writewav(spkr_signals,fs,fullpath_multiCh,writewavOpts);

if ~isempty(original)
    audiowrite([path ...
        orig_name '_' 'Original' ...
        orig_ext], original, fs, 'BitsPerSample',64);
end

end

