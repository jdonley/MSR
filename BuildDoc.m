%% Build documentation
tic;

wrkdir = 'M:\MSR\';
docdir = 'doc';

mainFile = '+Current_Systems\Evaluate_Current_System.m';
mainRunFile = { ...
    '+Current_Systems\loadCurrentSRsystem.m' ...
%     '+Results\Plot_Results.m' ...
};


docFiles = Tools.buildDocumentation( wrkdir, docdir, mainFile, {}, false, true );
docFiles = Tools.buildDocumentation( wrkdir, docdir, mainRunFile{1}, docFiles, true, true );
% docFiles = Tools.buildDocumentation( wrkdir, docdir, mainRunFile{2}, docFiles, true, true );


system('commit_doc.bat');
toc;