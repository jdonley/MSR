function [Path, err, room_details_path, RIR_Name__Details] = getRIRDatabasePath( setup, room, database_workingdir, method )
% Returns the path of the room impulse response database for a given setup and room or SR system
% 
% Syntax:	[Path, err, room_details_path, RIR_Name__Details] = ...
%               getRIRDatabasePath( ...
%                   setup, room, database_workingdir, method )
% 
% Inputs: 
% 	setup - Speaker_Setup.loudspeaker_setup object to generate RIR path from
% 	room - Room_Acoustics.Room object to generate RIR path from
% 	database_workingdir - The directory (folder) for the databases
% 	method - the method used to generate the path (for backwards
%            compatibility) (using old methods may result in ambiguities)
%            (defaults to the newest method)
% 
% Outputs: 
% 	Path - The room impulse response database path
% 	err - Error flag is false unless error is thrown when trying to
% 	      generate path
% 	room_details_path - A subpath generated from the room parameters only
% 	RIR_Name__Details - Part of the name of the final database file
% 
% Example: 
%     SYS = Current_Systems.loadCurrentSRsystem;
%     Room_Acoustics.getRIRDatabasePath( ...
%         SYS.Main_Setup(1), ...
%         SYS.Room_Setup, ...
%         SYS.system_info.Drive)
% 
% See also: loadCurrentSRsystem, SR_System, loudspeaker_setup, Room 

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015-2017
% Date: 13 November 2015
% Version: 0.1 (13 November 2015)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Settings

latest_method = 'new2';
if nargin < 4
    method = latest_method;
end
if nargin < 3
    database_workingdir = Settings.DefaultDatabasePath;
end


RIR_Database_Path = [database_workingdir ...
            '+Room_Acoustics\' ...
            '+RIR_Database\'];

err = false;
try
    
    if strcmpi(method, latest_method)
        
        [~,~,~,array_style_dir, spkr_array_dir, zone_positions_dir] = Soundfield_Database.getDatabasePath(setup, '', database_workingdir, 'new4');
        
        room_details_path = ['+' room.Room_Size_txt 'Dim_' ...
                             num2str(room.Wall_Absorb_Coeff) 'Ab\'];
        
        RIR_Name__Details = [num2str(room.NoReceivers) 'Rec_' ...
                             room.Reproduction_Centre_txt 'Ctr'];

        Path = [RIR_Database_Path ...
                array_style_dir, ...
                spkr_array_dir, ...
                zone_positions_dir, ...
                room_details_path, ...
                'RIRs__' RIR_Name__Details '.mat'];
    
    elseif strcmpi(method, 'new')
        
        SpeakerLayoutFolder = ['+' num2str(setup.Radius*2) 'm_SpkrDia_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];

        RIR_Database_Path = [RIR_Database_Path ...
            SpeakerLayoutFolder ...
            num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X,10)) 'Bx_' ...
            num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y,10)) 'By_' ...
            num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X,10))  'Qx_' ...
            num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y,10))  'Qy\'];
        
        RIR_Name__Details = [num2str(setup.Loudspeaker_Count) 'Src_' ...
            num2str(room.NoReceivers) 'Rec_' ...
            room.Room_Size_txt 'Dim_' ...
            room.Reproduction_Centre_txt 'Ctr_' ...
            num2str(room.Wall_Absorb_Coeff) 'Ab'];
        
        Path = [RIR_Database_Path 'RIRs__' RIR_Name__Details '.mat'];
        
    elseif strcmpi(method, 'old')
                
        RIR_Name__Details = [num2str(setup.Loudspeaker_Count) 'Src_' ...
            num2str(room.NoReceivers) 'Rec_' ...
            room.Room_Size_txt 'Dim_' ...
            room.Reproduction_Centre_txt 'Ctr_' ...
            num2str(room.Wall_Absorb_Coeff) 'Ab'];
        
        Path = [RIR_Database_Path 'RIRs__' RIR_Name__Details '.mat'];
            
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

