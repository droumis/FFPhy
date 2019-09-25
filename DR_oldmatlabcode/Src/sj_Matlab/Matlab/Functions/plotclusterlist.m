function [h] = plotclusterlist(aname, basedir, c, geom)
% [h] = plotclusterlist(adir, aname, basedir, celllist, geometry)
%       Plot, in a single figure, clusters for all of the place fields listed
%       in celllist. 
%
%	aname is the animal name (e.g. 'fred');
%
%	basedir should be the base directory where the animal's cluster data
%	are kept (e.g. '/dataxx/yourname/')
%
%	celllist is an n x 4 list of cells.  
%
%	geometry is [m n] and specifies the number of rows and columns for the
%	plot
%

% generate the subplot positions
rstart = 0.1;
rspace = .065;
cstart = 0.1;
cspace =  .15;
splotgeom = ones(geom(1), 1) * geom(2);
spos = getsubplotpos(splotgeom, 'rstart', rstart, 'rspace', rspace, ...
			'cstart', cstart, 'cspace', cspace);
a = animaldef(aname);

ncells = size(c,1);
spikes = loaddatastruct(a{2}, a{3}, 'spikes', c(1,1));

% go through the list of cells, exclude the exluded data and plot the place
% fields
i = 1;
for row = 1:geom(1)
    for col = 1:geom(2)
	if (i <= size(c,1))
	    subplot('Position', spos{row}{col});
	    tme = spikes{c(i,1)}{c(i,2)}{c(i,3)}{c(i,4)}.timerange / 10000;
	    [h1 h2] = plotcluster(aname, basedir, [2:5], tme, c(i,:));
	    f = title(sprintf('%d', i));
	    set(f, 'FontSize', 12);
	    p = get(f, 'Position');
            p(2) = p(2) * .925;
%            p(2) = p(2) * .9;
	    set(f, 'Position', p);
	    i = i + 1;
	end
    end
end
         
