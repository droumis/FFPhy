function out = getlindist(animaldir, animalprefix, epochs)
% out = getlindist(animaldir,animalprefix,epochs)
%
% Produces a cell structure with the fields:
% time, lindist
%   EPOCHS - N by 2 matrix, columns are [day epoch]
%

loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);

for i = 1:size(epochs,1)
    lindist = (linpos{epochs(i,1)}{epochs(i,2)}.statematrix.lindist);
    out{epochs(i,1)}{epochs(i,2)}.lindist= lindist;
    out{epochs(i,1)}{epochs(i,2)}.time = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;
end