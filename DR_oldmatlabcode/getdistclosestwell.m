function out = getdistclosestwell(animaldir, animalprefix, epochs)
%
%
%
loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);

for i = 1:size(epochs,1)
    [y ind ] = min(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.linearDistanceToWells,[],2);
    out{epochs(i,1)}{epochs(i,2)}.distwell= y;
    out{epochs(i,1)}{epochs(i,2)}.closestwell = ind;
    out{epochs(i,1)}{epochs(i,2)}.time = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;
end