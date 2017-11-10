function [Receiver_Signals] = Convolve_SpkrSigs_and_RIRs( Speaker_Signals, RIRs, Method )
% Convolves Room Impulse Responses (RIRs) with Loudspeaker Signals
%
% Syntax:	[Receiver_Signals] = CONVOLVE_SPKRSIGS_AND_RIRS( Speaker_Signals, RIRs )
%
% Inputs:
% 	Speaker_Signals - A 2D matrix of Speaker Signals such that each
%                     row represents a single loudspeaker signal.
%                     Speaker_Signals( NLoudspeakers, Signal_Length )
% 	RIRs - A 3D matrix of RIRs for each loudspeaker and receiver position
%          such that the first dimension is the receivers, the second is
%          the length of the RIR and the third is the loudspeakers.
%          RIRs( NReceivers, RIR_Length, NLoudspeakers )
%
% Outputs:
% 	Receiver_Signals - A 2D matrix comprising of the received signals
% 	defined by the length of the loudspeaker signals and the number of
% 	receivers deduced from the sizes of the input matrices.
%   Receiver_Signals( NLoudspeakers, NReceivers, Signal_Length + RIR_Length - 1 )
%
% Example:
%     SYS = Current_Systems.loadCurrentSRsystem;
%     SYS.signal_info.method = SYS.signal_info.methods_list{1};
%     spkrPath = Broadband_Tools.getLoudspeakerSignalPath( ...
%         SYS.Main_Setup(1), ...
%         SYS.signal_info, ...
%         SYS.system_info.LUT_resolution, ...
%         SYS.system_info.Drive);
%     spkrFiles = Tools.getAllFiles(spkrPath);
%     spkrSigs = audioread(spkrFiles{1});
%     DBpath = Soundfield_Database.getDatabasePath( ...
%         SYS.Main_Setup(1), ...
%         SYS.system_info.LUT_resolution, ...
%         SYS.system_info.Drive);
%     DB = load(DBpath);
%     receivedSigs = Convolve_SpkrSigs_and_RIRs( ...
%         spkrSigs, ...
%         DB.RIRs.Bright_RIRs, ...
%         'FFT');
%
% See also: getLoudspeakerSignalPath, getDatabasePath, loadCurrentSRsystem

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 04 August 2015
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default Method
if nargin < 3
    Method = 'builtin';
end

% Check matrix compatibility
if size(Speaker_Signals, 1) ~= size(RIRs, 3)
    error('Number of loudspeakers along the corresponding dimension does not match!')
end

NLoudspeakers = size(Speaker_Signals, 1);
NReceivers = size(RIRs,1);
Signal_Length = size(Speaker_Signals, 2);
RIR_Length = size(RIRs, 2);

Receiver_Signals_Length = Signal_Length + RIR_Length - 1;

if ~strcmp(Method,'FFT')
    Receiver_Signals = zeros(NLoudspeakers, NReceivers, Receiver_Signals_Length, 'like', Speaker_Signals);
    
    for spkr = 1:NLoudspeakers
        for rec = 1:NReceivers
            
            spkr_sig = Speaker_Signals(spkr,:);
            rir = RIRs(rec,:,spkr);
            
            if strcmp(Method,'builtin')
                Receiver_Signals( spkr, rec, : ) = conv(spkr_sig,rir);
            elseif strcmp(Method,'FFT_slow')
                Receiver_Signals( spkr, rec, : ) = Tools.fconv(spkr_sig.',rir.');
            end
            
        end
    end
else
    Receiver_Signals = zeros(NReceivers, Receiver_Signals_Length, 'like', Speaker_Signals);
    %parfor rec = 1:NReceivers
    for rec = 1:NReceivers
        NRec = length(rec);
        spkr_sig = repmat(Speaker_Signals.', 1, NRec);
        rir = reshape(permute(RIRs(rec,:,:),[3 1 2]), NRec * NLoudspeakers, RIR_Length).';
        
        rec_sig = Tools.fconv( spkr_sig, rir );
        %rec_sig = Tools.fconv( gpuArray(spkr_sig), gpuArray(rir) );
        
        Receiver_Signals(rec,:) = sum(rec_sig,2).';%reshape( rec_sig', NLoudspeakers, NRec, Receiver_Signals_Length);
    
    end
end

end
