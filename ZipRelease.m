%% Create release ZIP with Documentation
tic;

wrkdir = 'M:\MSR\';
docdir = 'doc';
zipFname = 'release';

mainFile = '+Current_Systems\Evaluate_Current_System.m';
mainRunFile = '+Current_Systems\loadCurrentSRsystem.m';


docFiles = Tools.buildDocumentation( wrkdir, docdir, mainFile, {}, false, true );
Tools.buildDocumentation( wrkdir, docdir, mainRunFile, docFiles, true, true );

AllDocFiles = Tools.getAllFiles([wrkdir docdir]);

relFiles = Tools.buildReleaseZIP( wrkdir, zipFname, mainFile, AllDocFiles, false );
Tools.buildReleaseZIP( wrkdir, zipFname, mainRunFile, relFiles, true );


toc;