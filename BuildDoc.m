%% Build documentation

wrkdir = 'M:\MSR\';
docdir = 'doc';

mainFile = '+Current_Systems\Evaluate_Current_System.m';
docFiles = Tools.buildDocumentation( wrkdir, docdir, mainFile, {}, false, true );

mainFile = '+Current_Systems\loadCurrentSRsystem.m';
Tools.buildDocumentation( wrkdir, docdir, mainFile, docFiles, true, true );
