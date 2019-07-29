function out = getpostlastwell(animaldir,animalprefix,epochs,varargin)
% out = getpostlastwell(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, postlast
% EPOCHS - N by 2 matrix, columns are [day epoch]
% intended to be used to exclude period after the last well visit
% post well may be contaminated with sleep periods or noisy periods when the animal
% is being transported to the sleep box

% Author Demetris Roumis June 2019

plotfigs = 0;
if ~isempty(varargin)
    assign(varargin{:});
end

loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);

for i = 1:size(epochs,1)
    timevec = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;
    out{epochs(i,1)}{epochs(i,2)}.time = timevec;
   
    % exlude times after last well visit
    welllastidx = find(diff( ...
        linpos{epochs(i,1)}{epochs(i,2)}.statematrix.wellExitEnter(:,2)),1,'last');
    [~,midx] = max(timevec);
    lastidx = min([midx welllastidx]);

    postlast = (1:length(timevec) > lastidx)';
    postlast_list = vec2list(postlast, timevec);

    out{epochs(i,1)}{epochs(i,2)}.postlast = postlast;
    out{epochs(i,1)}{epochs(i,2)}.postlast_list = postlast_list;
                
    if plotfigs
        figure
        pf = plot(postlast*100, 'y+', 'DisplayName','exl_prefirst');
        pf.MarkerSize = 8;
        axis tight
        pause
    end
end
