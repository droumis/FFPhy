% simclust = getsimclust(aname, basedir, linfields, taskstruct, clustcorrthresh, propthresh, placecorrthresh)
%	Returns a list of similar clusters for this animal based on 
%	1.  The correlation of 4D histograms of spike amplitudes across days
%	2.  The proportional between clusters.
%	3.  The correlation between placefields of putatively similar clusters 
%
%	aname is the name of the main animal directory (e.g. 'alex')
%
%	basedir is the name of the directory where that animal directory is
% 	located (e.g. '/data10/yourname/')
%
%	linfields is the linear place field structure containing all days of
%	data
%	
%	taskstruct is a single task structure with data from all days
%
%	clustcorrthresh is the minimum correlation between histograms that 
%	indicates possible correspondence. (suggested value = .25) 
%
%	propthresh is the proportion multipler for calculating the maximum 
%	distances between clusters on each dimension. Currently the code
%	expects two values [pt1 and pt2] where pt1 is the multipler for cases
%	where only 1 cluster is close and pt2 is for cases where > 1 cluster
%	is close.  Two clusters are considered close when the distance between
%	them is less that the pt1 or pt2 times the distance from the farthest
%	cluster to the origin. Suggested values: [0.1 0.2];
%	
%	placecorrthresh is a final test for the minimum correlation between 
%	the place fields of corresponding clusters.  If any pair of place 
%	fields are correlated at or above this level, the two clusters are 
%	considered to be corresponding. 
%	
function [simclust] = getsimclust(aname, basedir, linfields, task, clustcorrthresh, propthresh, placecorrthresh)

runs = [2 4 6];

simclust.aname = aname;
simclust.basedir = basedir;
simclust.clustcorrthresh = clustcorrthresh;
simclust.propthresh = propthresh;
simclust.placecorrthresh = placecorrthresh;

pushd(basedir);

filelist = {};
datasetnum = [];
tetnum = [];
nfiles = 0;


cd(aname);

% get a list of all of the dataset namesNow, this could be done easily with some nested for loops, but that really does violate the spirit in which such challenges are issued. I did not know of a function in MATLAB that did that, but I figured someone had invented that wheel before, so off to the File Exchange. A few minutes later I remembered the difference between combinations and permutations and had found two great new picks for this week.
dtmp = dir();
datadirs = dtmp(3:end);

% go into each data directory and each tetrode name and check for a matclust
% file.   

for d = 1:length(datadirs)
    cd(datadirs(d).name);
    dtmp = dir();
    tetdir = dtmp(3:end);
    for t = 1:length(tetdir)
	if (tetdir(t).isdir)
	    cd(tetdir(t).name);
	    mdir = dir();
	    % search for a file starting with matclust
	    for m = 3:length(mdir)
		if (strncmp(mdir(m).name, 'matclust', 8))
		    % there is a matclust file, so add it to the list
		    nfiles = nfiles + 1;
		    filelist{nfiles} = strcat(pwd, '/', mdir(m).name);
		    datasetnum(nfiles) = str2num(datadirs(d).name(end-1:end));
		    tetnum(nfiles) = str2num(tetdir(t).name(1:2));
		end
	    end
	    cd ..
	end
    end
    cd ..
end

% go through all combinations of files associated with the same tetrode 
% to get a set of correlation values
ind = 1;
simclust.corresp = struct([]);
for d1 = 1:nfiles
    % check to see if the task structure exists for this dataset and if not
    % move on to the next iteration of the loop
    try 
	tmptask = task{datasetnum(d1)};
    catch
	continue;
    end
    load(filelist{d1});
    clust{1} = clustattrib.clusters;
    parms{1} = clustdata.params;
    nclust(1) = length(clust{1});
    for d2 = (d1+1):nfiles
	% check to see if the task structure exists for this dataset and if not
	% move on to the next iteration of the loop
	try 
	    tmptask = task{datasetnum(d2)};
	catch
	    continue;
	end
	if (tetnum(d1) == tetnum(d2))
	    load(filelist{d2});
	    clust{2} = clustattrib.clusters;
	    parms{2} = clustdata.params;
	    nclust(2) = length(clust{2});
	    if (min(nclust) > 1)
		[c clustvect withindist clustdist] = clustsim(clust, parms);
		if (c > clustcorrthresh) 
		    % get the list of cluster vectors
		    cv1 = clustvect{1};
		    cv2 = clustvect{2};

		    % go through all pairs of clusters 
		    closeclust = zeros([2 nclust]);
		    totaldist = zeros([2 nclust]);
		    for c1 = 1:nclust(1)
			for c2 = 1:nclust(2)
			    % start with the proportion threshold for 1 
			    % close clusters 
			    for p = 1:2
				% get the individual distances
				d = clustdist(c1, c2, :);
				d = abs(d(:)');
				maxdist = max([cv1(c1,:) ; cv2(c2,:)]) * ...
				    propthresh(p); 
				% check to see if the distance along each dimension is less that
				% the cluster vector * prop
				if (sum(d < maxdist) == length(d))
				    closeclust(p, c1, c2) = 1;
				    totaldist(p, c1,c2) = sqrt(sum(d.^2));
				end
			    end
			end
		    end
		    % 
		    % go through and check to see if there are multiple clusters from one
		    % dataset that are close to a cluster from the other dataset, and if
		    % so, pick the cluster with the closest within cluster vector as the
		    % corresponding cluster
		    % start with a the rows
		    w1 = withindist{1};
		    w2 = withindist{2};
		    for p = 1:2
			tmpcc = closeclust(p,:,:);
			tmpcc = reshape(tmpcc, size(tmpcc,2), size(tmpcc,3)); 
			for r = 1:size(tmpcc, 1)
			    ctmp = find(tmpcc(r,:));
			    % get the distances between within cluster vectors
			    if (length(ctmp) > 1)
				% there are two or more close clusters in this row, so we find
				% the one that has the closest within cluster vector
				mindist = 100;
				for c = ctmp
				    tmpdist = dist(w1(r), w2(c));
				    if (tmpdist < mindist)
					mindist = tmpdist;
					mincol = c;
				    end
				end
				tmpcc(r,:) = 0;
				tmpcc(r,mincol) = 1;
				closeclust(p,r,:) = tmpcc(r,:);
			    end
			end
			% now do the columns
			for c = 1:size(tmpcc,2)
			    rtmp = find(tmpcc(:,c))';
			    % get the distances between within cluster vectors
			    if (length(rtmp) > 1)
				% there are two or more close clusters in this column, so we 
				% find the one that has the closest within cluster vector
				mindist = 100;
				for r = rtmp
				    tmpdist = dist(w1(r), w2(c));
				    if (tmpdist < mindist)
					mindist = tmpdist;
					minrow = r;
				    end
				end
				tmpcc(:,c) = 0;
				tmpcc(minrow,c) = 1;
				closeclust(p,:,c) = tmpcc(:,c);
			    end
			end
		    end
		    % now figure out which proportion threshold to use
		    if (length(find(closeclust(2,:,:))) > 1)
			% there are 2 or more close clusters using the higher
			% threshold, so we use that one.
			pthresh = 2;
		    elseif (length(find(closeclust(1,:,:))) == 1)
			pthresh = 1;
		    else
			pthresh = 0;
		    end
		    if (pthresh) 
			% initialize the same clust structure
			simclust.corresp(ind).dataset = ...
				    [datasetnum(d1) datasetnum(d2)];
			simclust.corresp(ind).tetnum = tetnum(d1);
			simclust.corresp(ind).cells = [];
			simclust.corresp(ind).corrruns = [];
			% get a list of environments
			tracklist = {};
			ttmp = task{datasetnum(d1)};
			for r = 1:length(ttmp)
                            if isfield(ttmp{r}, 'environment')
				tracklist{end+1} = ttmp{r}.environment;
			    end
			end
			ttmp = task{datasetnum(d2)};
			for r = 1:length(ttmp)
                            if isfield(ttmp{r}, 'environment')
				tracklist{end+1} = ttmp{r}.environment;
			    end
			end
			tracklist = unique(tracklist);
			% for each environment, find the epochs from each
			% dataset for that environment
			for t = 1:length(tracklist)
			    environ(t).d1 = [];
			    ttmp = task{datasetnum(d1)};
			    for r = 1:length(ttmp)
				if (isfield(ttmp{r}, 'environment') & ...
				   (strcmp(ttmp{r}.environment, tracklist{t})))
				    environ(t).d1(end+1) = r;
				end
			    end
			    environ(t).d2 = [];
			    ttmp = task{datasetnum(d2)};
			    for r = 1:length(ttmp)
				if (isfield(ttmp{r}, 'environment') & ...
				   (strcmp(ttmp{r}.environment, tracklist{t})))
				    environ(t).d2(end+1) = r;
				end
			    end
			end
			% get the place fields for each pair of two clusters and
			% compute the total correlation 
			tmpcc = closeclust(p,:,:);
			tmpcc = reshape(tmpcc, size(tmpcc,2), size(tmpcc,3)); 
			[c1, c2] = find(tmpcc);
			for c = 1:length(c1)
			    % go though all environments
			    for e = 1:length(environ)
				if (~isempty(environ(e).d1) & ...
				    ~isempty(environ(e).d2))
				    % go through all combinations of elements
				    r1 = environ(e).d1;
				    r2 = environ(e).d2;
				    [x y] = meshgrid(r1, r2);
				    x = x(:);
				    y = y(:);
				    for r = 1:length(x)
					try 
					    lf1 = linfields{datasetnum(d1)}{x(r)}{tetnum(d1)}{c1(c)};
					catch
					    break;
					end
					try 
					    lf2 = linfields{datasetnum(d2)}{y(r)}{tetnum(d2)}{c2(c)};
					catch
					    break;
					end
					if (~isempty(lf1) & ~isempty(lf2))
					    totalsame = 0;
					    len1 = min(length(lf1{1}(:,5)), length(lf2{1}(:,5)));
					    len2 = min(length(lf1{3}(:,5)), length(lf2{3}(:,5)));
					    tmp1 = [];
					    tmp2 = [];
					    for j = 1:4
						if (j <= 2)
						    l = len1;
						else
						    l = len2;
						end
						% get rid of NaNs
						lf1{j}(find(~isfinite(lf1{j}(1:l,5))), 5) = 0;
						lf2{j}(find(~isfinite(lf2{j}(1:l,5))), 5) = 0;
						tmp1 = [tmp1 ; lf1{j}(1:l,5)];
						tmp2 = [tmp2 ; lf2{j}(1:l,5)];
					    end
					    cc = corrcoef(tmp1, tmp2);
					    if (cc(1,2) > placecorrthresh) 
						% this pair is corresponding
						simclust.corresp(ind).cells =...
						    [simclust.corresp(ind).cells ; ...
							[c1(c) c2(c)]];
						simclust.corresp(ind).corrruns = ...
						    [simclust.corresp(ind).corrruns ; ...
						    cc(1,2) x(r) y(r)];
					    end
					end
				    end
				end
			    end
			end
			if (~isempty(simclust.corresp(ind).cells))
			    % we found one or more corresponding cells, so we
			    % save the totaldist matrix and move on to the next
			    % index 
			    tmptd = totaldist(p,:,:);
			    tmptd = reshape(tmptd,size(tmptd,2),size(tmptd,3)); 
			    simclust.corresp(ind).totaldist = tmptd;
                            ind = ind + 1;
			end
		    end
		end
	    end
	end
    end
end

% get rid of the final element if the cell list is empty
if ((length(simclust.corresp) > 0) & isempty(simclust.corresp(end).cells))
    simclust.corresp = simclust.corresp(1:(end-1));
end

popd
