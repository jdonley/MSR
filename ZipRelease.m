%% Create release ZIP with Documentation
tic;

wrkdir = 'M:\MSR\';
docdir = 'doc';
zipFname = 'release\release';

mainFile = '+Current_Systems\Evaluate_Current_System.m';
mainRunFile = {'+Current_Systems\loadCurrentSRsystem.m' ...
               '+Results\Plot_Results.m'};


docFiles = Tools.buildDocumentation( wrkdir, docdir, mainFile, {}, false, true );
docFiles = Tools.buildDocumentation( wrkdir, docdir, mainRunFile{1}, docFiles, true, true );
docFiles = Tools.buildDocumentation( wrkdir, docdir, mainRunFile{2}, docFiles, true, true );

AllDocFiles = Tools.getAllFiles([wrkdir docdir]);
AuxFiles = {'LICENSE'; 'README.md'};

relFiles = Tools.buildReleaseZIP( wrkdir, zipFname, mainFile, [AllDocFiles; AuxFiles], false );
relFiles = Tools.buildReleaseZIP( wrkdir, zipFname, mainRunFile{1}, relFiles, true );
relFiles = Tools.buildReleaseZIP( wrkdir, zipFname, mainRunFile{2}, relFiles, true );

cellfun(@(f) copyfile( f, 'release\mfiles\'), docFiles);

toc;