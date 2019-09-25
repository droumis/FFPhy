% [h1 h2] = plotclust(aname, basedir, parmcol, times, index)
%	Loads the matclust file for the cell specified by index, finds the two
%	dimensions with the largest values for the cell and plots the
%	cluster for the specified epoch in those dimensions using the times in 
%	the spikes structure.  
%
%	aname is the name of the main animal directory (e.g. 'alex')
%
%	basedir is the name of the directory where that animal directory is
% 	located (e.g. '/data10/yourname/')
%
%	parmcol is the columns for the parameters you want to search over
%	(usually 2-5 for amplitudes)
%
%	times is [starttime endtime] in seconds
%
%	index is [day epoch tetrode cell]

function [h1 h2] = plotcluster(aname, basedir, parmcol,tme, ind);

pushd(basedir);

filelist = {};
datasetnum = [];
tetnum = [];
nfiles = 0;
day = ind(1);
epoch = ind(2);

% replace the clist epoch with the epoch from the training filter

dirname = sprintf('%s/%s%02d', aname, aname, ind(1));
cd(dirname);

% get the tetrode directory
tetnumstr = sprintf('%02d*', ind(3));
dtmp = dir(tetnumstr);
cd(dtmp.name);
mdir = dir();
for m = 3:length(mdir)
    if (strncmp(mdir(m).name, 'matclust', 8))
	% there is a matclust file, so we load it up
	load(mdir(m).name);
    end
end
clust = clustattrib.clusters;
parms = clustdata.params;

% convert the times to timestamp units
tme = tme * 10000;

% get the events that are between the start and end times
validpind = find((parms(:,1) >= tme(1)) & (parms(:,1) <= tme(2)));
validp = parms(validpind, parmcol);

validc = intersect(clust{ind(4)}.index, validpind);


% find the two largest dimensions
cellparm = parms(validc, parmcol);
[y i] = sort(mean(cellparm));
dim = i(3:4);

% plot all of the points
mx = max(validp(:,dim));
h1 = plot(validp(:,dim(1)), validp(:,dim(2)), '.', 'Color', [.6 .6 .6]);
set(h1, 'MarkerSize', 1);
hold on

% plot the cluster
h2 = plot(cellparm(:,dim(1)), cellparm(:,dim(2)), 'b.');
set(h2, 'MarkerSize', 4);

a = axis;
set(gca, 'FontSize', 10, 'XLim', [0 mx(1)*1.05], 'YLim', [0 mx(2)*1.05]);

f = xlabel(sprintf('Amp. ch. %d (uV)', dim(1)));
set(f, 'FontSize', 10);
% move the label closer
p = get(f, 'Position');   
%p(2) = -mx(2) / 5;
%p(2) = -mx(2) / 8;
p(2) = -mx(2) / 6.5;
set(f, 'Position', p);

f = ylabel(sprintf('Amp. ch. %d (uV)', dim(2)));
set(f, 'FontSize', 10);
% move the label closer
%p = get(f, 'Position');   
%p(2) = p(2) * .9;
%set(f, 'Position', p);

%f = title(sprintf('cell %d %d %d %d', ind));
%p = get(f, 'Position');   
%p(2) = p(2) * .9;
%set(f, 'Position', p);


popd;
