function out = getsplittimes(animaldir, animalprefix, epochs, ntimeseg, timeseg)
% out = getalllinstate(animaldir,animalprefix,epochs, ntimeseg, timeseg)
% Produces a cell structure with the fields:
% time and includeseg
% EPOCHS - N by 2 matrix, columns are [day epoch]
%
% ntimeseg specifies that the number of equal sized time segments
% timeseg is the list of segment numbers that should be included
%

loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);
for i = 1:size(epochs,1)
    out{epochs(i,1)}{epochs(i,2)}.time = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;

    % get the set of times corresponding to each equally sized segment
    segtimes = linspace(out{epochs(i,1)}{epochs(i,2)}.time(1), ...  
			out{epochs(i,1)}{epochs(i,2)}.time(end), ntimeseg+1); 
    segind = lookup(segtimes, out{epochs(i,1)}{epochs(i,2)}.time);

    out{epochs(i,1)}{epochs(i,2)}.includeseg = ...
	zeros(size(out{epochs(i,1)}{epochs(i,2)}.time));

    % for each segment, put ones in the indeces corresponding to that segment
    % if it is listed in timeseg
    for s = 1:ntimeseg
	if (~isempty(find(timeseg == s)))
	    out{epochs(i,1)}{epochs(i,2)}.includeseg(segind(s):segind(s+1)) = 1;
	end
    end
end

