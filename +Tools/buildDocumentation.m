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

PubOpts ={...
    'format',               'html',             ...
    'stylesheet',           '',                 ...
    'createThumbnail',      true,               ...
    'figureSnapMethod',     'entireGUIWindow',  ...
    'imageFormat',          'png',              ...
    'maxHeight',            [],                 ...
    'maxWidth',             [],                 ...
    'useNewFigure',         true,               ...
    'evalCode',             false,              ...
    'catchError',           true,               ...
    ...'codeToEvaluate',       [],                 ...
    'maxOutputLines',       Inf,                ...
    'showCode',             false,              ...
    };

%% Get all dependencies of the main file
flist = matlab.codetools.requiredFilesAndProducts( MainFile );                         % Full dependency file list

%% Only keep files in working directory
docFiles = flist(contains(flist,WorkingDir));                                          % Files to document

%% Create documentation directory
docdirs = strrep(cellfun(@fileparts, strrep(docFiles, WorkingDir, ''),'un',0),'+',''); % Determine structure
htmldir = [WorkingDir DocDir filesep 'html' filesep];                                  % Set html directory
w = warning('query','MATLAB:MKDIR:DirectoryExists');                                   % Get mkdir warning state
warning('off','MATLAB:MKDIR:DirectoryExists');                                         % Turn off unnecessary warning
cellfun(@(a) mkdir( [htmldir a] ), unique(docdirs));                                   % Create html folder structure
warning(w.state,'MATLAB:MKDIR:DirectoryExists');                                       % Reset warning state

%% Find indices of compatible files only
[~,~,e]=cellfun(@fileparts, docFiles,'un',0);                                          % Get file extensions
Icode = cellfun(@(x) (strcmp(x,'.m') || strcmp(x,'.mlx')),e);                          % Determine if file is only code (determine if compatible for documentation)

%% Generate documentation
cellfun(@(f,d) publish(f,PubOpts{:},'outputDir',[htmldir d]),docFiles(Icode),docdirs(Icode),'un',0);

end
