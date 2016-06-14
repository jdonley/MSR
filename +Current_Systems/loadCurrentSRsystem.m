function SYS = loadCurrentSRsystem
%LOADCURRENTSYSTEM Loads the current Soundfield Reproduction system

System_Path = 'M:\MSR\+Current_Systems\'; 

Current_System = 'IEEETransactions_System_A';













FileExt = '.m';
SYS = loadFromName( [System_Path Current_System FileExt] );

end


function SYS = loadFromName( name )
[p,n,e] = fileparts(name);
fnDir = getClassDirs([p filesep]);
if fnDir ~= -1
    Sys = str2func([fnDir n]);
    SYS = Sys();
else
    fin = fopen(name);
    data = fread(fin); fclose(fin);
    tmpnm = [tempname e];
    [p_tmp,n_tmp] = fileparts(tmpnm);
    data = replaceFuncName(data,n_tmp);% alter function line
    fout = fopen(tmpnm,'w');
    fwrite(fout,data); fclose(fout);
    addpath(p_tmp);
    Sys = str2func(n_tmp);
    SYS = Sys();
    rmpath(p_tmp);
    delete(tmpnm);
end
end

function classDirs = getClassDirs(FullPath)
classDirs = '';
classes = strfind(FullPath,'+');
if ~isempty(classes)
    for c = 1:length(classes)
        clas = FullPath(classes(c):end);
        stp = strfind(clas,filesep);
        classDirs = [classDirs  clas(2:stp(1)-1) '.'];
    end
else
    classDirs = -1;
end
end

function data = replaceFuncName(data,newname)
NL{1} = char([13, 10]);
NL{2} = sprintf('\n');
indCRLF = strfind(data', NL{1});
if isempty(indCRLF)
    indCRLF = strfind(data', NL{2});
end
indFunc = strfind(lower(data'),'function');
indFuncEnd = indCRLF(indCRLF>indFunc);
funcLine = indFunc : indFuncEnd(1)-1;
firstLine = char(data(funcLine)');

ind_eq = indFunc + strfind(firstLine,'=');
if isempty(ind_eq)
    ind_eq = indFunc + length('function');
end
ind_br = indFunc + strfind(firstLine,'(');
if isempty(ind_br)
    ind_br = indFuncEnd;
end

data = [data(indFunc:ind_eq-1)' newname data(ind_br-1:end)']';
end




