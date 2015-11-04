function [DB, err] = loadRIRDatabaseFromSetup( setup, room, database_workingdir, method )
%GETDATABASEFROMSETUP Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    method = 'new';
end
if nargin < 3
    database_workingdir = 'Z:\';
end

RIR_Name__Details = [num2str(setup.Loudspeaker_Count) 'Src_' ...
    num2str(room.NoReceivers) 'Rec_' ...
    room.Room_Size_txt 'Dim_' ...
    room.Reproduction_Centre_txt 'Ctr_' ...
    num2str(room.Wall_Absorb_Coeff) 'Ab'];

RIR_Database_Path = [database_workingdir ...
            '+Room_Acoustics\' ...
            '+RIR_Database\'];

err = false;
try
    
    if strcmpi(method, 'new')
        
        SpeakerLayoutFolder = ['+' num2str(setup.Radius*2) 'm_SpkrDia_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];

        RIR_Database_Path = [RIR_Database_Path ...
            SpeakerLayoutFolder ...
            num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X,10)) 'Bx_' ...
            num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y,10)) 'By_' ...
            num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X,10))  'Qx_' ...
            num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y,10))  'Qy\'];
        
        DB = load([RIR_Database_Path 'RIRs__' RIR_Name__Details '.mat']);
        
    elseif strcmpi(method, 'old')
                
        DB = load([RIR_Database_Path 'RIRs__' RIR_Name__Details '.mat']);
            
    else
        error('Method to load RIR database from setup and room is not supported.')
    end
    
catch ex
    switch ex.identifier
        case 'MATLAB:load:couldNotReadFile'
            warning(['Could not load RIR database using the ' method ' method.']);
            err = true;
        otherwise
            rethrow(ex)
    end
end


end

