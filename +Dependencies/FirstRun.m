function err = FirstRun
%FIRSTRUN Checks for the existence of dependency folders and provides download URLs
%if they do not exist. If all dependency folders exist this function
%rewrites the Global_Settings file so FirstRun is false.

if Dependencies.CheckDependencies
    Dependencies.setFirstRunFalse
    err = false;
else
    err = true;
end

end

