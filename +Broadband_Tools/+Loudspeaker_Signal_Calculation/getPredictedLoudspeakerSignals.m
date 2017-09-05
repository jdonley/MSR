function [Loudspeaker_Signals, Original] = getPredictedLoudspeakerSignals( SYS )
%GETPREDICTEDLOUDSPEAKERSIGNALS Summary of this function goes here
%
% Syntax:	[Loudspeaker_Signals, Original] = ...
%               getPredictedLoudspeakerSignals( SYS )
%
% Inputs:
% 	SYS - Soundfield Reproduction system object
%
% Outputs:
% 	Loudspeaker_Signals - The loudspeaker signals predicted to cancel the
%                         soundfield from the received microphone signals
% 	Original - The orignal audio file
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
% Date: 31 August 2017
% Version: 0.1 (31 August 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy the system for dealing with the reception and transmission
mSYS = SYS; %  Microphone System (M-SYS)
lSYS = SYS; % Loudspeaker System (L-SYS)

%% Read Microphone Signals
mSYS.Main_Setup = mSYS.Main_Setup( ...
    strcmpi(SYS.signal_info.methods_list,'clean') );
mSYS.Room_Setup = mSYS.Room_Setup( ...
    strcmpi({mSYS.Room_Setup.SystemType},'receive') );

MicSigPath = Broadband_Tools.getMicrophoneSignalPath( mSYS );
MicSigFiles = Tools.keepFilesFromFolder( Tools.getAllFiles(MicSigPath), mSYS.signal_info.speech_filepath);

Original = audioread( ... 
    MicSigFiles{contains(lower(MicSigFiles),'original')});

MicSigFiles(~contains(MicSigFiles,'.mat'))=[];
MicSigFiles(~contains(MicSigFiles,mSYS.signal_info.input_filename))=[];

MSF = load(MicSigFiles{:}); % MicSigFile (MSF)

MicSigs = MSF.mic_signals;
Fs = MSF.fs;

%% Determine prediction length (hop size)
lSYS.Main_Setup = lSYS.Main_Setup( ...
    ~strcmpi(SYS.signal_info.methods_list,'clean')  );
lSYS.Room_Setup = lSYS.Room_Setup( ...
    strcmpi({lSYS.Room_Setup.SystemType},'transmit') );

Q = mSYS.Room_Setup.NoReceivers;
N = SYS.signal_info.Nfft;
c = SYS.signal_info.c;
buffLen = SYS.signal_info.predict_buff;
armethod = SYS.signal_info.AR_method;

spkLocsPol = lSYS.Main_Setup.Loudspeaker_Locations;    % Polar
micLocs = mSYS.Room_Setup.ReceiverPositions ...        % Cartesian
            - mSYS.Room_Setup.Reproduction_Centre([2 1 3]);
[spkLocs(:,1),spkLocs(:,2)] = ...
    pol2cart(spkLocsPol(:,1),spkLocsPol(:,2));

spkLocCent = mean(spkLocs,1);
micLocCent = mean(micLocs,1);

d = sum(abs(spkLocCent(1:2) - micLocCent(1:2)).^2).^.5;
hop = round(d/c*Fs);
ol = (N-hop)/N;

%%
if ol == 1
    Loudspeaker_Signals = MicSigs;
    return;
end

mics = 1:Q;
for mic_ = 1:numel(mics)
    mic = mics(mic_);
    x = MicSigs(:,mic);
    
    b = Tools.frame_data( [zeros( N*(buffLen-1),1);x], 1-(1-ol)/buffLen ,  N*buffLen);
    s = Tools.frame_data( [x; zeros( N*(buffLen-1),1)], ol ,  N);
    
    [~,sigPredicted(:,mic_)] = Broadband_Tools.PredictiveFraming( ...
        s,b, ...
        int64((1-ol)* N),...
        armethod);
    
end

Loudspeaker_Signals = sigPredicted(1:size(MicSigs,1),:);


end
