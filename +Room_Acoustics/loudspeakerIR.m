function [ ImpRsp, DirLst ] = loudspeakerIR( loudspeaker_setup, room, n_positions, n_samples, Fs, rec_positions, norm_positions )
%LOUDSPEAKERIR Computes impulse response from loudspeaker directivity
%according to the corresponding loudspeaker object.
%
%   This allows a single impulse response to be calculated so that
%   MCRoomSim can use the frequency dependent directivity along with the
%   image method to compute the room response of a directional source.
S = loudspeaker_setup;
if S.Loudspeaker_Count ~= 1
    error('Can only compute Impulse Response for a single loudspeaker in the given loudspeaker setup.')
end

freqs = linspace(0, Fs/2, n_samples+1);
freqs(end)=[];

xx = [];
yy = [];

if nargin >= 6 && ~isempty(norm_positions)
    xx = [norm_positions(:,1)] ...
          - room.Reproduction_Centre(1);
    yy = [norm_positions(:,2)] ...
          - room.Reproduction_Centre(2);
end

if nargin >= 5 && ~isempty(rec_positions)
    xx = [xx; ...
        [rec_positions.Bright_Receiver_Positions(:,1); ...
        rec_positions.Quiet_Receiver_Positions(:,1)] ...
        - room.Reproduction_Centre(1)];
    yy = [yy; ...
        [rec_positions.Bright_Receiver_Positions(:,2); ...
        rec_positions.Quiet_Receiver_Positions(:,2)] ...
        - room.Reproduction_Centre(2)];
end

% Some impule response relative positions
r = ones(1,n_positions);
th = linspace(-pi,pi, n_positions +1); th(end)=[];
[x,y] = pol2cart( th(:), r(:));

[spkrX, spkrY] = pol2cart( ...
    S.Loudspeaker_Locations(:,1), ...
    S.Loudspeaker_Locations(:,2));

xx = [xx, x + spkrX];
yy = [yy, y + spkrY];
zz = zeros(size(xx));

[~,H] = S.computeField( xx, yy, freqs );


% Normalise bright zone response if given some bright zone positions
if nargin >= 6 && ~isempty(norm_positions)
    H_EQ = mean( abs(H(1:size(norm_positions,1),:)) , 1); %Find normalise value from the positions given
    H = H(size(norm_positions,1)+1:end,:);  %Remove normalisation positions
    H(:,2:end) = H(:,2:end) ./ repmat(H_EQ(2:end),size(H,1),1); %Normalise
    H = H ./ (max(abs(H(:))) * loudspeaker_setup.Loudspeaker_Object.P1);
end

%Do IR stuff here
Z = cat(2, H, conj(flip(H,2)) );
[az, el, rr] = cart2sph( xx - spkrX, yy - spkrY, zz );
ir = real( ifft( Z, [], 2 ));
ir(isnan(ir))=0;
ir = flip( circshift(ir, size(ir,2)/2, 2 ), 2 );

ImpRsp = ir(:,end/2+1:end);
DirLst = [az(:)/pi*180, el(:)/pi*180, rr(:)];
end

