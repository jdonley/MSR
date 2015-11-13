function [DB, err] = loadDatabaseFromSetup( setup, database_res, database_workingdir, method )
%GETDATABASEFROMSETUP Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    method = 'new';
end

if nargin < 3
    database_workingdir = 'Z:\';
end

[path, err] = Soundfield_Database.getDatabasePath( setup, database_res, database_workingdir, method );

DB = load(path);    

end

