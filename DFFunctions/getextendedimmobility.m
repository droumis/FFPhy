function out = getextendedimmobility(animaldir,animalprefix,epochs,varargin)
% out = getextendedimmobility(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, immobile
% EPOCHS - N by 2 matrix, columns are [day epoch]
% intended to be used to exclude periods of extended immobility, based on
% distance traveled within a moving window. This allows for keeping
% immobile states that don't last long, i.e. at a reward well.. but exclude
% periods when the animal falls asleep. Used by DR to keep awake 'engaged'
% rips on the wtrack, but exclude times when he was in a different,
% 'disengaged' state.

% Author Demetris Roumis June 2019

pre_excl_win = 60; %seconds
post_excl_win = 10; %seconds
ax_srate = 30; %Hz
cm_thresh = 2; %cm
sig = 900;
gausnumpoints = 1200;
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
        linddist_win = conv(diff(lindist),ones(pre_excl_win*ax_srate,1),'same');
        
        linddist_win = zeros(length(lindist),1);
        linddist_win(2:end) = movsum(diff(lindist), ...
            [pre_excl_win*ax_srate post_excl_win*ax_srate],'omitnan');
        
        g = gaussian(sig,gausnumpoints);
        u = smoothvect(abs(linddist_win), g);
        exl_immobile =  u <= cm_thresh;
        excl_immobile_list = vec2list(exl_immobile, timevec);
        fprintf('%d %d :: %d intervals \n',epochs(i,1),epochs(i,2), length(excl_immobile_list));
        
        if plotfigs
            figure
            p = plot(exl_immobile*100, '+', 'DisplayName','exl_immobile');
            p.Color = [1 0 1 .01];
            p.MarkerSize = 10;
            hold on
            sm = plot(u, 'k');
            ld = plot(lindist, 'b');
            hold off
            axis tight
            pause
        end
        
        out{epochs(i,1)}{epochs(i,2)}.immobile = exl_immobile;
        out{epochs(i,1)}{epochs(i,2)}.immobile_list = excl_immobile_list;
    catch
        fprintf('skipping %s linpos day %d epoch %d \n', animaldir, epochs(i,1), ...
            epochs(i,2));
        continue
    end
end
