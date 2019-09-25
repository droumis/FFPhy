function out = get2dstate(animaldir,animalprefix,epochs)
% out = get2dstate(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, headdir, velocity
% EPOCHS - N by 2 matrix, columns are [day epoch]
%

loaddays = unique(epochs(:,1));
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
for i = 1:size(epochs,1)
    if size(pos{epochs(i,1)}{epochs(i,2)}.data,2) > 5
        tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,8);
    else
        tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,5);
    end
    timestep = (pos{epochs(i,1)}{epochs(i,2)}.data(2,1) - pos{epochs(i,1)}{epochs(i,2)}.data(1,1));
    out{epochs(i,1)}{epochs(i,2)}.time = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
    out{epochs(i,1)}{epochs(i,2)}.headdir = pos{epochs(i,1)}{epochs(i,2)}.data(:,4);
    out{epochs(i,1)}{epochs(i,2)}.velocity = tmpvelocity;
    out{epochs(i,1)}{epochs(i,2)}.immobilitytime = cumsum_reset(tmpvelocity<2)*timestep;
    
end
%--------------------------------------------------------------------------

function out = cumsum_reset(inputvect)

resetpoints = find(diff(inputvect) < 0)+1;
out = cumsum(inputvect);
for i = 1:length(resetpoints)
    out(resetpoints(i):end) = out(resetpoints(i):end)-out(resetpoints(i)-1);
end