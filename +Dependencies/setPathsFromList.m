function setPathsFromList
%SETPATHSFROMLIST Sets and saves the dependency paths

userpath('M:\MSR');
addpath(Dependencies.DependencyList);
rmpath(Dependencies.FolderExceptionList);
savepath;

end

