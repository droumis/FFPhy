function out = getpriortofirstwell(animaldir,animalprefix,epochs,varargin)
% out = getpriortofirstwell(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, immobile
% EPOCHS - N by 2 matrix, columns are [day epoch]
% intended to be used to exclude period before the first well visit
% often this is a very noisy period that contains the time when the animal
% is being transported from the sleep box to the well, and therefore has
% invalid position data

% Author Demetris Roumis June 2019

plotfigs = 0;

if ~isempty(varargin)
    assign(varargin{:});
end

loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);

for i = 1:size(epochs,1)
    try
        timevec = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;
        out{epochs(i,1)}{epochs(i,2)}.time = timevec;
        
        lindist = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.linearDistanceToWells(:,1);

        % exlude times before first well visit
        wellfirstidx = find(diff( ...
            linpos{epochs(i,1)}{epochs(i,2)}.statematrix.wellExitEnter(:,2)),1);
        centerwell_enteridx = find(lindist < 4, 1);
        outerwell_enteridx = find(lindist > max(lindist)-4 , 1);
        firstidx = min([outerwell_enteridx centerwell_enteridx wellfirstidx]);
        
        excl_prefirst = (1:length(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time) ...
            < firstidx)';
        excl_prefirst_list = vec2list(excl_prefirst, ...
            linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time);
        
        
        fprintf('%d %d :: %d intervals \n',epochs(i,1),epochs(i,2), ...
            length(excl_prefirst_list));
        
        if plotfigs
            figure
            pf = plot(excl_prefirst*100, 'y+', 'DisplayName','exl_prefirst');
            pf.MarkerSize = 8;
            axis tight
            pause
        end
        
        out{epochs(i,1)}{epochs(i,2)}.prefirst = excl_prefirst;
        out{epochs(i,1)}{epochs(i,2)}.prefirst_list = excl_prefirst_list;
    catch
        fprintf('skipping %s linpos day %d epoch %d \n', animaldir, epochs(i,1), ...
            epochs(i,2));
        continue
    end
end
