function [ AllDependenciesExist ] = CheckDependencies
%CHECKDEPENDENCIES Checks for the dependency folders
AllDependenciesExist = false;

global Settings
Dpath = Settings.DependenciesPath;
sep = ';';

d = Dependencies.Dependencies;

ReqDepFolders = cellfun(@(x) [Dpath x sep],...
    {d(find(~cellfun(@isempty,{d.Required}))).Name},...
    'UniformOutput',false);

ReqDepURLs = {d(find(~cellfun(@isempty,{d.Required}))).URL};

ReqDepNames = {d(find(~cellfun(@isempty,{d.Required}))).Name};

isMissing = false;
for f = 1:numel(ReqDepFolders)
    if ~exist(strrep(ReqDepFolders{f},';',''),'dir')
        isMissing = true;
        Tools.simpleWarning({...
            ['The dependency "' ReqDepNames{f} ...
            '" was not found in the folder "' strrep(strrep(ReqDepFolders{f},';',''),'\','\\') '"'], ...
            'Try obtaining the dependency from: ' });
        Tools.printHyperlink( ReqDepURLs{f} );
        fprintf('\n');
    end
end

if isMissing
    Tools.simpleWarning({...
            ['Please install all missing dependencies and restart MATLAB.\n\n'] });
else
    AllDependenciesExist = true;
end

end

