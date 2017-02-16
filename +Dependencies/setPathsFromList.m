function setPathsFromList
%SETPATHSFROMLIST Sets and saves the dependency paths
global Settings

if (isfield(Settings, 'FirstRun') && Settings.FirstRun) || ...
        ~isfield(Settings, 'FirstRun')
   isError = Dependencies.FirstRun;
   if isError
       return;
   end
end
if Settings.ChangeUserpath
    userpath(Settings.ToolboxPath);
end
addpath(Dependencies.DependencyList);
rmpath(Dependencies.FolderExceptionList);
savepath;

end

