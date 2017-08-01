function [DB, err, DB_path] = loadRIRDatabaseFromSetup( setup, room, database_workingdir, method )
%GETDATABASEFROMSETUP Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    method = 'new2';
end
if nargin < 3
    database_workingdir = 'Z:\';
end

[DB_path, err] = Room_Acoustics.getRIRDatabasePath(setup,room,database_workingdir, method);

DB = load(DB_path);

end

