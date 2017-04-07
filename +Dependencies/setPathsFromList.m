function setPathsFromList
%SETPATHSFROMLIST Sets and saves the dependency paths
global Settings
All_Settings.Global_Settings

if (isfield(Settings, 'FirstRun') && Settings.FirstRun) || ...
        ~isfield(Settings, 'FirstRun')
   isError = Dependencies.FirstRun;
   if isError
       return;
   end
end
if Settings.ChangeUserpath
    userpath(Settings.ToolboxPath);
    cd(strrep(userpath,';','')); % Required from R2017a onwards
end
addpath(Dependencies.DependencyList);
rmpath(Dependencies.FolderExceptionList);
savepath;

end

