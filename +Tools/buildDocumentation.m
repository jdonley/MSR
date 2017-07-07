function buildDocumentation( WorkingDir, DocDir, MainFile )
%BUILDDOCUMENTATION Generates documentation HTML and builds MATLAB search database for dependencies of a main file
% 
% Syntax:	BUILDDOCUMENTATION( WORKINGDIR, DOCDIR, MAINFILE )
% 
% Inputs: 
% 	WorkingDir  - The working directory of the project to document
% 	zipFileName - The directory for the documentation
% 	MainFile    - The path to the file for which to determine the dependencies
%                 and consequently build the documentation
% 
% Example: 
%     wrkdir = 'C:\myProject\';
%     docdir = 'doc';
%     mainFile = 'C:\myProject\examples\fullyfledgedexample.m';
%     buildDocumentation( wrkdir, docdir, mainFile )
% 
% See also: publish, builddocsearchdb, doc, requiredFilesAndProducts

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 07 July 2017 
% Version: 0.1 (07 July 2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get all dependencies of the main file
flist = matlab.codetools.requiredFilesAndProducts( MainFile );                         % Full dependency file list

%% Only keep files in working directory
docFiles = flist(contains(flist,WorkingDir));                                          % Files to document

%% Create documentation directory and html files
newclsdirs = unique(cellfun(@fileparts, strrep(docFiles, WorkingDir, ''),'un',0));     % Determine structure
newdirs = strrep(newclsdirs,'+','');                                                   % Rename class directories
htmldir = [WorkingDir DocDir filesep 'html' filesep];                                  % Set html directory
cellfun(@(a) mkdir( [htmldir a] ), newdirs);                                           % Create html folder structure



end
