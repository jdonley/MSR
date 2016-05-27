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
        audioSig = audioSig ./ rms(audioSig(:));
    catch err
        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
            continue; % Skip unsupported files
        end
    end
    speech = [speech; audioSig];
end

% spect_ = periodogram(speech,[],nfft,fs);
% spect_(2:end-1) = spect_(2:end-1) / 2;
% %spect_ = sqrt( spect_ * fs * nfft );
% spect_ = sqrt( spect_ * fs * length(speech) );


%%%%%%%%%%%%%%%%%%%%%%%%%%
win_=rectwin(nfft);
ovlap = 0;
speech_framed = buffer(speech,nfft);

% win_=hamming(nfft);
% ovlap = 0.5;
% speech_framed = enframe(speech,win_,nfft*ovlap,'z').';

x=speech;
[pxx,frqs_]=pwelch(x,win_,nfft*ovlap,nfft,fs,'power');


% xdft_Tot = zeros(1,nfft/2+1);
% for fr = 1:size(speech_framed,2)
%  x=speech_framed(:,fr);
% 
% U = sum(win_)^2;
% B = size(speech_framed,2);
% N = length(speech);
% N_B = nfft;
% K = length(frqs_);
% dK = mean(diff(frqs_));
% [n,k]=ndgrid(0:N_B-1,frqs_);
% x_ = repmat(x,1,K);
% xdft_(frqs_/dK + 1) = sum( x_ .* exp(  -2*1j*pi*k.*n / (N_B*dK)  ), 1 );
% xdft_ = abs(xdft_).^2;
% xdft_Tot = xdft_Tot + xdft_;
% end
% xdft_=xdft_Tot / (B * N_B^2);
% xdft_(2:end-1) = 2 * xdft_(2:end-1);
% xdft_=sqrt(xdft_);


% figure(1)
% % plot(frqs_,mag2db(spect_)); 
% plot(frqs_,mag2db(xdft_));hold on;
% plot(frqs_,pow2db(pxx));
% set(gca,'XScale','log');
% xlim([50 10000]);
% grid on
% hold off;

spect_=pxx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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