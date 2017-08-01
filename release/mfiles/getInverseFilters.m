function [ invIRs ] = getInverseFilters( IRs, invIRlen, f_band, fs, reg, method )
%GETINVERSEFILTERS 
if nargin < 6
    method = 'toeplitz';
end

Nspkrs = size(IRs,2);
invIRs = zeros(invIRlen*fs,Nspkrs);

for spkr = 1:Nspkrs
    %fprintf('\n\t processing loudspeaker %d ... ',spkr);
    if strcmpi(method, 'kirkeby')
        %invIRs(:,spkr) = Tools.IRcompactingKirkebyFilter(IRs(:,spkr), invIRlen, f_band, fs, reg);
        invIRs(:,spkr) = Tools.invFIR('complex', IRs(:,spkr), size(IRs,1), 0, invIRlen*fs, f_band, reg, 0);
    
    elseif strcmpi(method, 'toeplitz')
        invIRs(:,spkr) = Tools.invimplms( IRs(:,spkr), invIRlen*fs, 0);
    
    end
    %fprintf('completed');
end

%fprintf('\n');
end

