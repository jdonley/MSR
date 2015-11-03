function [DB, err] = loadDatabaseFromSetup( setup, database_res, database_workingdir, method )
%GETDATABASEFROMSETUP Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    method = 'new';
end

if nargin < 3
    database_workingdir = 'Z:\';
end



    err = false;
    try
        
        if strcmpi(method, 'new')
            
            Database_Path = [database_workingdir '+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
            DB = load([Database_Path 'Database_' ...
                num2str(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Angle/pi*180) 'degB_' ...
                num2str(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Angle/pi*180) 'degQ_' ...
                num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' ...
                database_res '.mat']);
            
        elseif strcmpi(method, 'old')
            
            Database_Path = [database_workingdir '+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
            DB = load([Database_Path 'LUT_Weight_vs_Frequency_' num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' database_res '.mat']);
            
        elseif strcmpi(method, 'old_zones_swapped')
            
            Database_Path = [database_workingdir '+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
            DB = load([Database_Path 'LUT_Weight_vs_Frequency_' num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' database_res '_zones_swapped.mat']);
            
        else
            error('Method to load database from setup is not supported.')
        end
        
    catch ex
        switch ex.identifier
            case 'MATLAB:load:couldNotReadFile'
                warning(['Could not load database using the ' method ' method.']);
                err = true;
            otherwise
                rethrow(ex)
        end
    end
    

end

