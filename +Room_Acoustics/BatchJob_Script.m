clear;close all;fclose all;clc;

% % ROOM 1
% % Anechoic
Room_Size = [10 10 10]; %Anechoic
Wall_Absorption_Coeff = 1.0;

% % ROOM 2
% % Small Open Plan Office
% %Room_Size = [5 6 4];  %Small Open Plan Office
% %Wall_Absorption_Coeff = 0.3;

% % ROOM 3
% % Medium Open Plan Office
% Room_Size = [8 10 3];   %Medium Open Plan Office
% Wall_Absorption_Coeff = 0.3;

% % ROOM 4
% % Cafe / Restaurant
%Room_Size = [9 14 3];   %Cafe/Restaurant
%Wall_Absorption_Coeff = 0.3;

% % ROOM 5
% % Small Open Plan Office
%Room_Size = [4 9 3];   %Small Open Plan Office
%Wall_Absorption_Coeff = 0.3;

%%
setups = [1 2 3 4 5];
for scheme = setups
    if scheme == 1
        % % Setup and Privacy Scheme 1
        mask_type = 'FlatMask';
        pw_angle = 0;
    elseif scheme == 2
        % % Setup and Privacy Scheme 2
        mask_type = 'FlatMask';
        pw_angle = 15;
    elseif scheme == 3
        % % Setup and Privacy Scheme 3
        mask_type = 'ZoneWeightMask';
        pw_angle = 15;
    elseif scheme == 4
        % % Setup and Privacy Scheme 4
        mask_type = 'FlatMask';
        pw_angle = 90;
    elseif scheme == 5
        % % Setup and Privacy Scheme 5
        mask_type = 'ZoneWeightMask';
        pw_angle = 90;
    end
    
    Room_Acoustics.Apply_RIRs.Reverberant_MSR_Analysis_batchfunc(Room_Size, Wall_Absorption_Coeff, mask_type, pw_angle);
    
end