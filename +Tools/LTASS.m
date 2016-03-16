function [ spect, frqs ] = LTASS( folder, nfft )
%LTASS Computes the Long-Term Average Speech Spectrum from a folder of
%speech files
% 
% Syntax:	[ spect ] = LTASS( folder )
% 
% Inputs: 
% 	folder - The path to the folder containing the speech files
% 
% Outputs: 
% 	spect - The LTASS spectrum
% 	frqs - The frequency vector for the spectrum
 
% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016
% Date: 2 March 2016
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Just incase this function tries to call getAllFiles within a class folder we
% should create a function handle for getAllFiles regardless
inf = dbstack('-completenames');
funcName = 'getAllFiles';
funcPath = inf.file;
classDirs = getClassDirs(funcPath);
getAllFiles_ = str2func([classDirs funcName]);

%% Start LTASS
if nargin < 2
    nfft = 1024;
end
files = getAllFiles_(folder);
speech=[];
F = length(files);
for file = 1:F
    try
        [audioSig,fs] = audioread(files{file});
    catch err
        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
            continue; % Skip unsupported files
        end
    end
    speech = [speech; audioSig];
end

spect_ = periodogram(speech,[],nfft,fs);
spect_(2:end-1) = spect_(2:end-1) / 2;
spect_ = sqrt( spect_ * fs * nfft );
% nfft = length(speech);
% spect_full = fft(speech);
% spect_ = abs( spect_full(1:nfft/2+1) );
frqs_ = (0:fs/nfft:fs/2)';

% frqs = [0, logspace(log10(frqs_(2)),log10(fs/2),nfft/2)]';
% spect = interp1( frqs_, spect_, frqs );
frqs = frqs_;
spect = spect_;

end

function classDirs = getClassDirs(FullPath)
    classDirs = '';
    classes = strfind(FullPath,'+');
    for c = 1:length(classes)
        clas = FullPath(classes(c):end);
        stp = strfind(clas,filesep);
       classDirs = [classDirs  clas(2:stp(1)-1) '.'];
    end
end