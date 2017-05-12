%% Create release ZIP

wrkdir = 'M:\MSR\';
zipFname = 'release';
mainFile = '\+Current_Systems\Evaluate_Current_System.m';

Tools.buildReleaseZIP( wrkdir, zipFname, mainFile )
