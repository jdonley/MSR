function deleteResultsFile( resultspath, result_types )
%DELETERESULTSFILE Summary of this function goes here

pExts = {'.csv','.mat'}; %possible file extensions

for i=1:length(result_types)
    filenm = [resultspath result_types{i} '_Results'];
    for p=1:numel(pExts)
        filenme = [filenm,pExts{p}];
        if exist(filenme,'file')
            delete(filenme);
            wrnCol = [255,100,0]/255;
            cprintf(wrnCol, ['The previous ' result_types{i} ' result file was deleted.\n'] );
        end
    end
end

end

