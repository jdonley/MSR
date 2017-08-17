function saveMicrophoneSignals( path, name, ext, mic_signals, mic_count, original, orig_name, orig_ext, fs )
% Summary of this function goes here
% 
% Syntax:	saveMicrophoneSignals( path, name, ext, mic_signals, mic_count, original, orig_name, orig_ext, fs )
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% Outputs: 
% 	output1 - Description
% 	output2 - Description
% 
% Example: 
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 17 August 2017
% Version: 0.1 (17 August 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(mic_count)
    mic_count = size(mic_signals,2);
end

if ~exist(path,'dir'); mkdir(path); end

fullpath_multiCh = [path name num2str(mic_count) 'Ch' ext];
audiowrite(fullpath_multiCh, mic_signals, fs, 'BitsPerSample',64);

if ~isempty(original)
    audiowrite([path ...
        orig_name '_' 'Original' ...
        orig_ext], original, fs, 'BitsPerSample',64);
end

end

