function buildReleaseZIP( WorkingDir, zipFileName, MainFile )
%BUILDRELEASEZIP Creates a ZIP file of all release dependencies for a main file
% 
% Syntax:	BUILDRELEASEZIP( WORKINGDIR, ZIPFILENAME, MAINFILE )
% 
% Inputs: 
% 	WorkingDir  - The working directory of the project to release
% 	zipFileName - The name the final ZIP file
% 	MainFile    - The path to the file for which to determine the dependencies
%                 and consequently build the release ZIP
% 
% Example: 
%     wrkdir = 'C:\myProject\';
%     zipFname = 'myProject_release';
%     mainFile = 'C:\myProject\examples\fullyfledgedexample.m';
%     buildReleaseZIP( wrkdir, zipFname, mainFile )
% 
% See also: requiredFilesAndProducts, zip

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 07 May 2017
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get all dependencies of the main file
flist = matlab.codetools.requiredFilesAndProducts( MainFile );                      % Full dependency file list

%% Only keep files in working directory
relFiles = flist(contains(flist,WorkingDir));                                       % Release files

%% Copy to temp directory so zipping retains folder structure
tmpdir = [tempname filesep];                                                        % Get temp directory
newdirs = unique(cellfun(@fileparts, strrep(relFiles, WorkingDir, ''),'un',0));     % Determine structure
cellfun(@(a) mkdir( [tmpdir a] ), newdirs);                                         % Create temp folder structure
cellfun(@(a,b) copyfile(a,[tmpdir b]), relFiles, strrep(relFiles, WorkingDir, '')); % Copy

%% Zip the entire folder
zip(zipFileName, tmpdir, WorkingDir);

%% Clean up
rmdir(tmpdir,'s');                                                                  % Delete all the temporarily copied files

end
