function out = sj_get2dstate(animaldir,animalprefix,epochs)
% out = get2dstate(animaldir,animalprefix,epochs)

% See kk_get2dstate also. Do a more accurate immobilitytime by using timestep 
% And use an absolure threshold for speed - should it be 3 or 5 for run,
% and 2 for sleep

% Produces a cell structure with the fields:
% time, headdir, velocity
% EPOCHS - N by 2 matrix, columns are [day epoch]
%

loaddays = unique(epochs(:,1));
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
for i = 1:size(epochs,1)
    timestep = (pos{epochs(i,1)}{epochs(i,2)}.data(2,1) - pos{epochs(i,1)}{epochs(i,2)}.data(1,1));
    posdata = pos{epochs(i,1)}{epochs(i,2)}.data;    
    if size(posdata,2)>5 % already smoothed position and filtered velocity
        tmpvelocity = abs(posdata(:,9)); % Should already be absolute value
        out{epochs(i,1)}{epochs(i,2)}.headdir = pos{epochs(i,1)}{epochs(i,2)}.data(:,8);

    else
        tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,5);
        out{epochs(i,1)}{epochs(i,2)}.headdir = pos{epochs(i,1)}{epochs(i,2)}.data(:,4);
    end
    % Using original position and velocity
    tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,5);
    out{epochs(i,1)}{epochs(i,2)}.headdir = pos{epochs(i,1)}{epochs(i,2)}.data(:,4);
    
    out{epochs(i,1)}{epochs(i,2)}.time = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
    out{epochs(i,1)}{epochs(i,2)}.velocity = tmpvelocity;
    %out{epochs(i,1)}{epochs(i,2)}.immobilitytime = cumsum_reset(tmpvelocity<(max(tmpvelocity)*.05))/30;
    %out{epochs(i,1)}{epochs(i,2)}.immobilitytime = cumsum_reset(tmpvelocity<(max(tmpvelocity)*.05))*timestep;
    out{epochs(i,1)}{epochs(i,2)}.immobilitytime = cumsum_reset(tmpvelocity<2)*timestep; % For sleep - speed 2
    
    
end
%--------------------------------------------------------------------------

function out = cumsum_reset(inputvect)

resetpoints = find(diff(inputvect) < 0)+1;
out = cumsum(inputvect);
for i = 1:length(resetpoints)
    out(resetpoints(i):end) = out(resetpoints(i):end)-out(resetpoints(i)-1);
end