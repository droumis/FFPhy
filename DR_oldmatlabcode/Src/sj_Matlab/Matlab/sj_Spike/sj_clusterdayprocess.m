function sj_clusterdayprocess(fpre, adir, apre, days, paramnum, varargin)
%clusterdayprocess(FOLDERPREFIX, ANIMALDIRECT, ANIMALPREFIX, DAYS, PARAMNUMS,
%                  OPTIONS)
%
% sj: Folder prefix is folder name with raw data. Get directory name for day by seraching 
% sj: 2-5 paramnum is amplitude, 14-17 is Wave1PC, 10-13 is Valley
%
%FOLDREPREFIX -- if the folders are called conley01, conley02, ... then 
%	FOLDERPREFIX is 'conley'
%ANIMALDIRECT -- the path to where the animal's processed data will be stored 
%	example '/data99/student/Con'
%ANIMALPREFIX -- also the first three letters of the animal's name (example 
%	'con') which will attach to the beginning of the .mat files containing 
%	the variables.
% DAYS -- the days to process
% PARAMNUMS -- the indeces of the matclust parameter numbers to be used for the
% 	cluster quality calculation (2-5 for the amplitude on the four channels)
%
% 	This function goes through all all tetrodes from the specified animal 
% `	and, for each epoch, calculates the L-ratio and IsolationDistance 
%	measures for every cluster.
%
%	The measures are saved in a single file ending in 'clustqual'
% 
%	Each dataset must have already been processed by 'dayprocess' to
%	produce the spike files which are used to delineate the epoch times


filelist = {};
datasetnum = [];
tetnum = [];
nfiles = 0;

if (fpre(end) ~= '/')
    fpre = [fpre,'/']; 
end

for i = days
    % load the spike data for this day
    spikes = loaddatastruct(adir, apre, 'spikes', i);
    if (isempty(spikes))
	continue;
    end
    % search through the spike struct for a non-empty structure for each epoch
    % this is a bit inefficient but it works
    nepochs = length(spikes{i});
    timerange = zeros(nepochs,2);
    for e = 1:nepochs
	for t = 1:length(spikes{i}{e})
	    for c = 1:length(spikes{i}{e}{t})
		if (~isempty(spikes{i}{e}{t}{c}))
		    timerange(e,:) = spikes{i}{e}{t}{c}.timerange;
		end
	    end
	end
    end

    if i < 10
       daystring = ['0',num2str(i)];
    else
       daystring = num2str(i);
    end
    
    % Get dirctory for current day
    %daydirect = [fpre,daystring];
    daydirname=[]; % initialize for each day to catch errors
    cd(fpre);
    dirtmp = dir();
    daydir = dirtmp(3:end);
    for dr = 1:length(daydir)
        currname = daydir(dr).name;
        if strcmp(currname(1:2),daystring)==1;
            daydirname = daydir(dr).name;
            break
        end
    end
    daydirect = [fpre,daydirname];
    % Done - directory for day
    disp(['Day ',daystring]);
    
    cd(daydirect);
    dtmp = dir();
    tetdir = dtmp(3:end);
    for t = 1:length(tetdir)
	if (tetdir(t).isdir)
	    cd(tetdir(t).name);
	    % get the first two numbers as the tetrode number
	    tetnum = str2num(tetdir(t).name(1:2));
	    mdir = dir();
	    % search for a file starting with matclust
	    for m = 3:length(mdir)
		if (strncmp(mdir(m).name, 'matclust', 8))
		    % there is a matclust file, so load it.
		    load(mdir(m).name);
		    clust = clustattrib.clusters;
		    nclust = length(clust);
		    % find any columns of the parameters that are all zeros and
		    % get rid of them
		    zcol = find(~sum(clustdata.params));
		    tmpparmnum = setdiff(paramnum, zcol);
		    parms = clustdata.params(:,tmpparmnum);
		    % go through the epochs and the clusters and compute 
		    % L-ratio and Isolation distance
		    for e = 1:nepochs
			eind = find((clustdata.params(:,1) >= ...
				timerange(e,1)) & ...
			            (clustdata.params(:,1) <= ...
				timerange(e,2)));
			for c = 1:nclust
			    % find the spike indeces that are in this epoch.
			    cind = intersect(eind, clust{c}.index);
			    if (length(cind) > 4)
				% subtract off the first eind so that we can use
				% cind to index into eind.
				cind = cind - eind(1) + 1;
				clustqual{i}{e}{tetnum}{c}.lratio = ...
				    L_Ratio(parms(eind,:), cind);
				clustqual{i}{e}{tetnum}{c}.isoldist = ...
				    IsolationDistance(parms(eind,:), cind);
			    else
				clustqual{i}{e}{tetnum}{c}.lratio = NaN;
				clustqual{i}{e}{tetnum}{c}.isoldist = NaN;
			    end
			end
		    end
		end
	    end
	    cd ..
	end
    end
    cd ..
    % save the data
    fname = sprintf('%sclustqual%02d', apre, i);
    fname = fullfile(adir, fname);
    save(fname, 'clustqual');
end

