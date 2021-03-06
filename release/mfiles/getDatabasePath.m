function [ Path, err, Sub_Path, array_style_dir, spkr_array_dir, zone_positions_dir, virtual_source_dir, database_filename ] = getDatabasePath( setup, database_res, database_workingdir, method )
%GETDATABASEPATH Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4 || isempty(method)
    method = 'new4';
end

if nargin < 3 || isempty(database_workingdir)
    database_workingdir = 'Z:\';
end

Sub_Path = '';
array_style_dir = '';
spkr_array_dir = '';
zone_positions_dir = '';
virtual_source_dir = '';

    err = false;
    try
                
%%%%%%%%%%%%%%%%%%%%%%%%% NEW4 METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
% Added virtual source parameters descriptions to dedicated directory
        if strcmpi(method, 'new4')
            sc = '_'; %Separation Character
            
            %%% Type of database
            pre_dir = ['+' 'Soundfield_Database' '\'];
            
            %%% Speaker array type, perpendicular distance to the centre
            %%% speaker, angle to the the centre speaker.
            array_style_dir = ['+' num2str(setup.Dimensionality) 'D' sc ...
                                   upper(setup.Speaker_Array_Type) 'array' sc ...
                                   num2str(setup.Radius) 'mPerpDist' sc ...
                                   num2str(setup.Speaker_Array_Centre) 'degCentre' '\'];
            
            %%% Number of loudspeakers and type of loudspeakers
            spkrTypeString = strrep(setup.Loudspeaker_Type, ' ','');
            spkr_type_dirstr = ['+' num2str(setup.Loudspeaker_Count) spkrTypeString 'Spkrs' sc];
            
            %%% Speaker array style.
            switch lower(setup.Speaker_Array_Type)
                case 'circle'
                    spkr_array_dir  = [spkr_type_dirstr num2str(setup.Speaker_Arc_Angle/180*pi * setup.Radius) 'mLen' '\'];
                    
                case 'line'
                    spkr_array_dir  = [spkr_type_dirstr ...
                        num2str( (setup.Loudspeaker_Count-1)  * (setup.Speaker_Spacing + setup.Loudspeaker_Dimensions(1)) ) 'mLen' '\'];
                    
                case '2line'
                    spkr_array_dir  = [spkr_type_dirstr ...
                        num2str( (ceil(setup.Loudspeaker_Count/2)-1)  * (setup.Speaker_Spacing + setup.Loudspeaker_Dimensions(1)) ) 'mLen' '\'];
                otherwise
                    error('Loudspeaker array type (circle, line, etc) from the given loudspeaker setup object does not have a database path defined.')
            end
            
            %%% Cartesian coordinates of bright zone and quiet zone
            zone_positions_dir = ['+' num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X*1e10)/1e10) 'Bx' sc ...
                                      num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y*1e10)/1e10) 'By' sc ...
                                      num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X*1e10)/1e10)  'Qx' sc ...
                                      num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y*1e10)/1e10)  'Qy' '\'];
                                  
            %%% Bright zone virtual source angle, distance and type
            virtual_source_dir = ['+' 'VSrc' sc ...
                                      setup.Multizone_Soundfield.Bright_Zone.SourceType sc ...
                                      num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg' sc ...
                                      num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Distance) 'm' '\' ];
                                              
            %%% Database resolution
            database_filename = ['Database' sc database_res '.mat'];
                              
            Sub_Path = [array_style_dir, ...
                    spkr_array_dir, ...
                    zone_positions_dir, ...
                    virtual_source_dir];
                
            Path = [database_workingdir, ...
                    pre_dir, ...
                    Sub_Path, ...
                    database_filename];  
                
%%%%%%%%%%%%%%%%%%%%%%%%% NEW3 METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
% Filepath now includes all parameters used within a loudspeaker setup
        elseif strcmpi(method, 'new3')
            
            %%% Type of database
            pre_dir = '+Soundfield_Database\';
            
            %%% Speaker array type, perpendicular distance to the centre
            %%% speaker, angle to the the centre speaker.
            array_style_dir = ['+' upper(setup.Speaker_Array_Type) 'array_' num2str(setup.Radius) 'mPerpDist_' num2str(setup.Speaker_Array_Centre) 'degCentre\'];
            
            %%% Number of loudspeakers and speaker array length.
            switch setup.Speaker_Array_Type
                case 'circle'
                    spkr_array_dir  = ['+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle/180*pi * setup.Radius) 'mLen\'];
                    
                case 'line'
                    spkr_array_dir  = ['+' num2str(setup.Loudspeaker_Count) 'Spkrs_' ...
                        num2str( (setup.Loudspeaker_Count-1)  * (setup.Speaker_Spacing + setup.Loudspeaker_Dimensions(1)) ) 'mLen\'];
                otherwise
                    error('Loudspeaker array type (circle, line, etc) from the given loudspeaker setup object does not have a database path defined.')
            end
            
            %%% Cartesian coordinates of bright zone and quiet zone
            zone_positions_dir = ['+' num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X,10)) 'Bx_' ...
                                      num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y,10)) 'By_' ...
                                      num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X,10))  'Qx_' ...
                                      num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y,10))  'Qy\'];
                                  
            %%% Bright zone virtual source angle, Database resolution
            database_filename = ['Database_' ...
                                  num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' ...
                                  database_res '.mat'];
                              
            Sub_Path = [array_style_dir, ...
                    spkr_array_dir, ...
                    zone_positions_dir];
                
            Path = [database_workingdir, ...
                    pre_dir, ...
                    Sub_Path, ...
                    database_filename];                  
                                  
%%%%%%%%%%%%%%%%%%%%%%%%% NEW2 METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
% Added support for line array database paths
% Improved zone position support in database path (full cartesian coords)        
        elseif strcmpi(method, 'new2')
            
            switch setup.Speaker_Array_Type
                case 'circle'
                    Database_Path = [database_workingdir ...
                        '+Soundfield_Database\' ...
                        '+' num2str(setup.Radius*2) 'm_ArrDia\' ...
                        '+' num2str(setup.Loudspeaker_Count) 'Spkrs_' ...
                        num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
                case 'line'
                    Database_Path = [database_workingdir ...
                        '+Soundfield_Database\' ...
                        '+' num2str( ...
                        (setup.Loudspeaker_Count-1)  * (setup.Speaker_Spacing + setup.Loudspeaker_Dimensions(1)) ...
                        ) 'm_ArrLen\' ...
                        '+' num2str(setup.Loudspeaker_Count) 'Spkrs_' ...
                        num2str(setup.Speaker_Spacing) 'Spac\'];
                otherwise
                    error('Loudspeaker array type does not have a database path defined.')
            end

            Path = [Database_Path ...
                '+' num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X,10)) 'Bx_' ...
                num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y,10)) 'By_' ...
                num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X,10))  'Qx_' ...
                num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y,10))  'Qy\' ...
                'Database_' ...
                num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' ...
                database_res '.mat'];
            
%%%%%%%%%%%%%%%%%%%%%%%%% NEW METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
% Added support for zone angle in database path
        elseif strcmpi(method, 'new')
                        
            Database_Path = [database_workingdir '+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
            Path = [Database_Path 'Database_' ...
                num2str(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Angle/pi*180) 'degB_' ...
                num2str(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Angle/pi*180) 'degQ_' ...
                num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' ...
                database_res '.mat'];
            
%%%%%%%%%%%%%%%%%%%%%%%%% OLD METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmpi(method, 'old')
            
            Database_Path = [database_workingdir '+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
            Path = [Database_Path 'LUT_Weight_vs_Frequency_' num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' database_res '.mat'];
            
        elseif strcmpi(method, 'old_zones_swapped')
            
            Database_Path = [database_workingdir '+Soundfield_Database\+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc\'];
            Path = [Database_Path 'LUT_Weight_vs_Frequency_' num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'deg_' database_res '_zones_swapped.mat'];
            
        else
            error('Method to get database path from setup is not supported.')
        end
        
    catch ex
        switch ex.identifier
            case 'MATLAB:load:couldNotReadFile'
                warning(['Could not get database path using the ''' method ''' method.']);
                err = true;
            otherwise
                rethrow(ex)
        end
    end

end

