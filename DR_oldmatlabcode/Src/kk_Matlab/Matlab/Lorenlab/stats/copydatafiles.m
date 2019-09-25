function copydatafiles(startdir, enddir)

cd(startdir)
dirnames = dir;
for i = 3:length(dirnames)
    if (dirnames(i).isdir && (length(dirnames(i).name) == 3))
        mkdir(enddir,dirnames(i).name);
        cd(dirnames(i).name);
        filenames = dir;
        for j = 3:length(filenames)
            if (~filenames(j).isdir)
          
                success = copyfile(filenames(j).name, fullfile(enddir,dirnames(i).name));
                disp([fullfile(dirnames(i).name,filenames(j).name),' ', num2str(success)]);
                
            end
        end
    end
    cd(startdir)
end