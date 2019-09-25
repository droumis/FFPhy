% pos = GETSUBPLOTPOS(nplotsperrow, options)
%       Given the nplotsperrow array, returns the coordinates for the subplots in a two
%       level cell where pos{1}{1} is in the upper left subplot
%       options are
%       'rspace', n    the spacing between rows (0 < n < 1)
%       'cspace', n    the spacing between columns (0 < n < 1)
%       'rstart', n    the offset of the rows relative to zero (0 < n < 1)
%       'cstart', n    the offset of the columns relative to zero (0 < n < 1)
function [pos] = getsubplotpos(n, varargin)

% assign the options
rspace = .05;
cspace = .05;
rstart = .1;
cstart = .1;
if (~isempty(varargin))
	assign(varargin{:});
end

nrows = length(n);
pos = cell(nrows,1);

rsize = (1 - (nrows - 1) * rspace - rstart - .075) / nrows;

% go through each of the rows and figure out the size of each plot
for i = 1:nrows
	ncols = n(i);
	csize = (1 - (ncols - 1) * cspace - cstart -.05) / ncols;
	pos{i} = cell(ncols,1);
	for j = 1:n(i)
		pos{i}{j} = [(cstart + (cspace+csize)*(j-1)) (rstart+(nrows-i)*(rsize+rspace)) csize rsize]; 
	end
end
