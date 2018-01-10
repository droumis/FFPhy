function out = excludenoiseevents(animaldir,animalprefix, epochs, eventconsname, tetfilter, varargin)
loaddays = unique(epochs(:,1));
eventcons = loaddatastruct(animaldir, animalprefix, eventconsname, loaddays);

for i = 1:size(epochs,1)
    
    %disp(sprintf('kk_getconstimes: (START) %s %d %d',animalprefix,epochs(i,1),epochs(i,2)))
    
    % empty error checking: possibly no eventcons..(perhaps not computed because of too much noise this epoch) if so, skip
    if isempty(eventcons)
        disp(sprintf('d %d e %d skipped, noise struct is empty',epochs(i,1),epochs(i,2)))
        out{epochs(i,1)}{epochs(i,2)}.time = nan;
        out{epochs(i,1)}{epochs(i,2)}.noise = nan;
        continue
    elseif (epochs(i,2) > length(eventcons{epochs(i,1)})) || isempty(eventcons{epochs(i,1)}{epochs(i,2)})
        disp(sprintf('d %d e %d skipped, no eventcons data',epochs(i,1),epochs(i,2)))
        out{epochs(i,1)}{epochs(i,2)}.time = nan;
        out{epochs(i,1)}{epochs(i,2)}.noise = nan;
        continue
    end
    
    % Apply tetfilter to find which entry (grouping of tetrodes) in eventcons is the right one.
    % If a number, then it's the index.
    % If a string, then search for the TF (tetfilter) index within eventcons
    if isnumeric(tetfilter)
        TF = tetfilter;
    else
       error('bum bum bummmmmmm')
    end
    
         
    
    ec = eventcons{epochs(i,1)}{epochs(i,2)}{TF};
    times = ec.timerange(1):1/ec.samprate:ec.timerange(end);
    ectimes = [ec.starttime(:) ec.endtime(:)];
    noise = list2vec(ectimes,times)';
    noisetime = (sum(noise)/1500/60);
    totaltime = (ec.timerange(end) - ec.timerange(1))/60;
    disp(sprintf('day%d ep%d -- %.02f of %.02f minutes labeled as noise', epochs(i,1),epochs(i,2), noisetime, totaltime))
    out{epochs(i,1)}{epochs(i,2)} = ec;
    out{epochs(i,1)}{epochs(i,2)}.time = times;
    out{epochs(i,1)}{epochs(i,2)}.noise = noise;
end

end