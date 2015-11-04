function deleteResultsFile( resultspath, info_filename, result_types )
%DELETERESULTSFILE Summary of this function goes here
%   Detailed explanation goes here

Num_ = strfind(info_filename,'__')+2;
Num_p = Num_(3);
z_t = sscanf([' ' info_filename(Num_p:end)],'%*[^_]_%[A-Za-z]');
ind = strfind(info_filename, z_t);
[ pw_a, ~, wei, m_t ] = Results.getInfoFromFilename(info_filename(1:ind-2));

for i=1:length(result_types)
    csvFilename = ([resultspath result_types{i} '_Results_' num2str(pw_a) 'pwAngle_' num2str(wei) 'weight' m_t '.csv']);    
    if exist(csvFilename,'file')    
        delete(csvFilename);
        warning( ['The previous ' result_types{i} ' result file was deleted.'] );
    end
end

end

