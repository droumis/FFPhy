function renameeegfilevariable(directoryname,fileprefix,days)

daytetlist = [];

for day = days
    % create the list of files for this day that we should filter
       tmpflist = dir(sprintf('%s/EEG/*thetagnd%02d-*.mat', directoryname, day));
       flist = cell(size(tmpflist));
       for i = 1:length(tmpflist)
	   flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
       end

    % go through each file in flist and filter it
    for fnum = 1:length(flist)
	% get the tetrode number and epoch
	% this is ugly, but it works
	dash = find(flist{fnum} == '-');
	epoch = str2num(flist{fnum}((dash(1)+1):(dash(2)-1)));
	tet = str2num(flist{fnum}((dash(2)+1):(dash(2)+3)));
	   
	%load the eeg file
	load(flist{fnum});
    
    %%MANUALLY CHANGE VARIABLE NAME HERE%%
    
    try
    thetagnd=theta;
    
    %% MANUALLY CHANGE BELOW AS WELL
    % save the resulting file
	thetafile = sprintf('%s/EEG/%sthetagnd%02d-%d-%02d.mat', ...
			    directoryname, fileprefix, day, epoch, tet);
        save(thetafile, 'thetagnd');
    disp('file saved!')
    catch
        disp('already correct')
    end
    end
end
