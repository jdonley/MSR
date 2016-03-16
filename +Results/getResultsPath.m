function [Path, err] = getResultsPath( setup, database_res, room, signal_info, database_workingdir, method )
%GETDATABASEFROMSETUP Summary of this function goes here
%   Detailed explanation goes here
latest_method = 'new';
if nargin < 6
    method = latest_method;
end
if nargin < 5
    database_workingdir = 'Z:\';
end
if nargin < 4
    signal_info = [];
end


Recordings_Path = [ ...
    database_workingdir ...
    '+Results\'];

err = false;
try
    
    if strcmpi(method, latest_method)
        sc = '_';
        
        [~,~,~,~,reproduction_info_dirs] = Broadband_Tools.getLoudspeakerSignalPath( setup, signal_info, database_res );
        
        spkr_sig_info_dir = [ ...
            '+' num2str(signal_info.f_low ) 'Hz-' ...
            num2str(signal_info.f_high) 'Hz' sc sc ...
            'method' sc signal_info.method filesep ];
        
        [~,~,room_info_dir1,room_info_dir2] = Room_Acoustics.getRIRDatabasePath( setup, room );
        room_info_dirs = [room_info_dir1, room_info_dir2, filesep];
        
        if isfield(signal_info, 'recording_type') && strcmpi(signal_info.recording_type, 'realworld')
            realworld_path = '+Physical_World\';
        else
            realworld_path = [];
        end
        
        Path = [Recordings_Path ...
            realworld_path, ...
            reproduction_info_dirs, ...
            room_info_dirs, ...
            spkr_sig_info_dir];
        
        
    elseif strcmpi(method, 'old')
        spkr_sig_dir = ['+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc_LUT_' database_res '\'];
                
        Path = [ ...
            Recordings_Path
            '+Reverb__' num2str(room.NoReceivers) 'Rec_' ...
            room.Room_Size_txt 'Dim_' ...
            room.Reproduction_Centre_txt 'Ctr_' ...
            num2str(room.Wall_Absorb_Coeff) 'Ab\' ...
            spkr_sig_dir];
            
    else
        error('Method to load Results path from setup and room is not supported.')
    end
    
catch ex
    switch ex.identifier
        case 'MATLAB:load:couldNotReadFile'
            warning(['Could not load Results path using the ' method ' method.']);
            err = true;
        otherwise
            rethrow(ex)
    end
end


end

