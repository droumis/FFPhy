function [h] = plofilterfields2d(f, fi, c, geom)
% [h] = plottrainingfields(trainingfilter, trainingindex, celllist, geometry)
%       Plot, in a single figure, heat maps for all of the place fields listed
%       in celllist. This uses the excludetimes from the training filter
%	The trainingindex should be [animal epoch1 epoch2] where
%		animal is the animal number for the trainingdata
%		epoch1 is the index of the epochfilter
%		epoch2 is the row corresponding to the day and epoch of the 
%			epochfilter
%	geometry is [m n] and specifies the number of rows and columns for the
%	plot
%


% get the day, animal and epoch
animal = f(fi(1)).animal;
ind = f(fi(1)).epochs{fi(2)}(fi(3),:)
day = ind(1);
epoch = ind(2);

% check to make sure that the day matchs
if ((c(1,1) ~= day))
    error('the day for the celllist must match the day and epoch for the training filter');
end

% generate the subplot positions
rstart = 1 / (40 * geom(1))
rspace = 2 * rstart;
cstart = 1 / (40 * geom(2))
cspace = 2 * cstart;
splotgeom = ones(geom(1), 1) * geom(2);
spos = getsubplotpos(splotgeom, 'rstart', rstart, 'rspace', rspace, ...
			'cstart', cstart, 'cspace', cspace);
% get the excludetime list
excludetime = f(fi(1)).excludetime{fi(2)}{fi(3)};

% load the spike and posifion data
spikes = loaddatastruct(animal{2}, animal{3}, 'spikes', day);
pos = loaddatastruct(animal{2}, animal{3}, 'pos', day);

ncells = size(c,1);

goodpos = ~isExcluded(pos{day}{epoch}.data(:,1), excludetime);
p = pos{day}{epoch}.data(goodpos,:);
% go through the list of cells, exclude the exluded data and plot the place
% fields
i = 1;
for row = 1:geom(1)
    for col = 1:geom(2)
	if (i <= size(c,1))
	    goodspikes = ~isExcluded(spikes{day}{epoch}{c(i,3)}{c(i,4)}.data(:,1), excludetime);
	    s = spikes{day}{epoch}{c(i,3)}{c(i,4)}.data(goodspikes,:);
	    rmap = openfieldrate(s, p, 2, 4);
	    
	    subplot('Position', spos{row}{col});
	    plotratemap(rmap.smoothedspikerate);
	    i = i + 1;
	end
    end
end
         
