function [DB, err] = loadRIRDatabaseFromSetup( setup, room, database_workingdir, method )
%GETDATABASEFROMSETUP Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    method = 'new';
end
if nargin < 3
    database_workingdir = 'Z:\';
end

[RIR_DB_path, err] = Room_Acoustics.getRIRDatabasePath(setup,room,database_workingdir, method);

DB = load(RIR_DB_path);

end

