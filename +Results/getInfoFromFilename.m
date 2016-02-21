function [ PW_angle, noise_mask, weight, mask_type ] = getInfoFromFilename( filename )
%GETINFOFROMFILENAME Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%
% % DEPRECATED  % %
%%%%%%%%%%%%%%%%%%%

% Work out Planewave Angle
pos = strfind(filename,'Hz_');
PW_angle = sscanf(filename(pos:end),'%*[^0-9^-]%d%*s');

if ~isempty(strfind(filename,'with'))
    % Work out noise mask
    pos = strfind(filename,'Angle');
    noise_mask = sscanf(filename(pos:end),'%*[^0-9^-]%d%*s');
    
    % Work out weight
    pos = strfind(filename,'dB');
    weight = sscanf(filename(pos:end),'%*[^0-9]%d%*s');
    
    % Work out mask type
    pos = strfind(filename,'weight');
    mask_type = sscanf(filename(pos:end),'weight%[^0-9]');
    mask_type = ['__' strrep(mask_type,'_','')];
else
    % Work out weight
    pos = strfind(filename,'Angle');
    weight = sscanf(filename(pos:end),'%*[^0-9]%d%*s');
    
    % Assign mask type and noise mask
    mask_type = 'NoMask';
    noise_mask = -Inf;
end

end

