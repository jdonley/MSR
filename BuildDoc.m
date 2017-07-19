%% Build documentation
tic;

wrkdir = 'M:\MSR\';
docdir = 'doc';

mainFile = '+Current_Systems\Evaluate_Current_System.m';
mainRunFile = '+Current_Systems\loadCurrentSRsystem.m';


docFiles = Tools.buildDocumentation( wrkdir, docdir, mainFile, {}, false, true );
Tools.buildDocumentation( wrkdir, docdir, mainRunFile, docFiles, true, true );


system('commit_doc.bat');
toc;