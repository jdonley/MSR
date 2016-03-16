function [Path, err] = getRealRecordingsPath( setup, database_res, room, signal_info, database_workingdir, method )
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


RealWorldRecordings_Path = '+Physical_World\';

err = false;
try
    
    if strcmpi(method, latest_method)
        
        [~,~,rec_pre_path,recording_subpath] = Results.getRecordingsPath( ...
            setup, ...
            database_res, ...
            room, ...
            signal_info, ...
            database_workingdir, ...
            method);
        
        Path = [rec_pre_path, RealWorldRecordings_Path, recording_subpath];
                   
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

