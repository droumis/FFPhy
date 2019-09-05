function out = excludenoiseevents(animaldir,animal, epochs, eventconsname, varargin)

tetfilter = 1;
exclpad = 0; % in seconds. padding pre and post event start/end to also exclude
stdthresh = 0;
excludeman = 0;
if ~isempty(varargin)
    assign(varargin{:})
end

if excludeman
    fprintf('excludeman\n');
    manNoiseEvents=load_data('filterframework','noiseEvents', animal, 'animpos', 0);
end

loaddays = unique(epochs(:,1));
eventcons = loaddatastruct(animaldir, animal, eventconsname, loaddays);

for i = 1:size(epochs,1)
    
    %disp(sprintf('kk_getconstimes: (START) %s %d %d',animalprefix,epochs(i,1),epochs(i,2)))
    
    % empty error checking: possibly no eventcons..(perhaps not computed because of too much noise this epoch) if so, skip
    if isempty(eventcons)
        fprintf('d %d e %d skipped, noise struct is empty\n',epochs(i,1),epochs(i,2))
        out{epochs(i,1)}{epochs(i,2)}.time = nan;
        out{epochs(i,1)}{epochs(i,2)}.noise = nan;
        continue
    elseif (epochs(i,2) > length(eventcons{epochs(i,1)})) || isempty(eventcons{epochs(i,1)}{epochs(i,2)})
        fprintf('d %d e %d skipped, no eventcons data\n',epochs(i,1),epochs(i,2))
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
    if stdthresh
        y = ec.maxthresh > stdthresh;
        fprintf('%d %d: %d noise events >%dstd \n', epochs(i,1), epochs(i,2), sum(y), ...
            stdthresh)
        ec.maxthresh = ec.maxthresh(y);
        ec.starttime = ec.starttime(y);
        ec.endtime = ec.endtime(y);
        ec.peak_value = ec.peak_value(y);
        ec.peak_index = ec.peak_index(y);
        ec.total_area = ec.total_area(y);
        ec.area_midpoint_index = ec.area_midpoint_index(y);
        ec.energy_midpoint_index = ec.energy_midpoint_index(y);
        ec.absolute_maxthresh = ec.absolute_maxthresh(y);
    end
    
    times = ec.timerange(1):1/ec.samprate:ec.timerange(end);
    if excludeman
        if ~isempty(manNoiseEvents.events) && any(ismember(manNoiseEvents.events(:,[1 2]), epochs(i,[1 2]), 'rows'))
            evidxEp = ismember(manNoiseEvents.events(:,[1 2]), epochs(i,[1 2]), 'rows');
            manEvStart = manNoiseEvents.events(evidxEp, 3);
            ec.starttime = unique([ec.starttime; manEvStart]);
            ec.endtime = unique([ec.endtime; manEvStart]);
        end
        ec.excludeman = manNoiseEvents;
    end
    ectimes = [ec.starttime(:)-exclpad ec.endtime(:)+exclpad];
    noise = list2vec(ectimes,times)';
    noisetime = (sum(noise)/ec.samprate/60);
    totaltime = (ec.timerange(end) - ec.timerange(1))/60;
    fprintf('day%d ep%d -- %.02f of %.02f minutes noise\n', epochs(i,1),epochs(i,2), ...
        noisetime, totaltime)
    out{epochs(i,1)}{epochs(i,2)} = ec;
    out{epochs(i,1)}{epochs(i,2)}.time = times;
    out{epochs(i,1)}{epochs(i,2)}.noise = noise;
    
end

end