function ingIntervals = getIngestionIntervals(pos, dio, task, day, epochs)
% get intervals when animal is ingesting reward
% @DKR-2021

% rew out was 1 sec after IN. so default -1.1 for 100 ms prior to IN.
rewOutOffset4Start = -1.1; % seconds. 
leavingVelThresh = 4;

ingIntervals = {};
for e = 1:length(epochs)
    posfields = strsplit(pos{day}{epochs(e)}.fields);
    posidx = getFieldIndex(posfields, {'time' 'vel-loess'});
    posTimeVel = pos{day}{epochs(e)}.data(:, posidx);
    
    % get times of reward dio output
    rewT = get_rewardTimes(task, dio, day, epochs(e));
    arrivs = [rewT(:,1)+rewOutOffset4Start rewT(:,2)];
    leaving = [];
    
    for i = 1:size(rewT, 1)
        posIdxAtRewT = lookup(rewT(i,1), posTimeVel(:,1));
        velAtWellPosIdx = find(posTimeVel(posIdxAtRewT:end,2) > ...
            leavingVelThresh, 1);
        movingAwayPosIdx = velAtWellPosIdx + posIdxAtRewT;
        try
            leaving = [leaving; posTimeVel(movingAwayPosIdx, 1)];
        catch
            fprintf('out of bounds\n');
            arrivs(i,:) = [];
            continue
        end
    end
    
    ingIntervals{e,1} = [arrivs(:,1) leaving];
end
end
