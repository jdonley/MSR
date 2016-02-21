function deleteResultsFile( resultspath, result_types )
%DELETERESULTSFILE Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(result_types)
    csvFilename = ([resultspath result_types{i} '_Results.csv']);    
    if exist(csvFilename,'file')    
        delete(csvFilename);
        warning( ['The previous ' result_types{i} ' result file was deleted.'] );
    end
end

end

