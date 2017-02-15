Time_Lapse_Path = 'D:\Pictures\Speaker_Array_Timelapse\';
files=Tools.getAllFiles(Time_Lapse_Path);
p=[];
imDate={};
for f=1:length(files)
    p(f) = ~isempty(strfind(files{f},'.JPG'));
    if p(f)
    imDate{f}=imfinfo(files{f});
    imDate{f}=imDate{f}.DateTime;
    end
end
files = files(logical(p));
[~,I]=sort(imDate(logical(p)));
files = files(I);
F = length(files);
for f=1:F
   
    [fpath,fname,fext]=fileparts(files{f});
    movefile( ...
        files{f},...
        [Time_Lapse_Path, 'IMG_', num2str(f), fext]);
end