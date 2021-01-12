function rewTimes = get_rewardTimes(task, DIO, day, epochs)


outdios = task{day}{epochs(1)}.outputdio;
rewTimes = [];
for e = 1:length(epochs)
    rewTimes = [rewTimes; get_dios_times(outdios, DIO, day, epochs(e))];
end

% get rid of duplicate timestamps in output.. (trodes bug)
% [R, uid, ~] = unique(rewTimes(:,1));
% R = [R rewTimes(uid,2)];
% [~, sidx] = sort(R(:,1));
% rewTimes = R(sidx,:);
% rewTimes = rewTimes(2:end,:); % drop artificial first

% confirm that the animal is at the reward well for each rewTime
% condition: animal pos is within X of well at time of rew
% condition: there is a well input trigger 1 second preceding rew time
% rewTimes = rewTimes(distanceFromWellAtRewTime < 20);
% posAtRewTime = pos(lookup(rewTimes(:,1), pos(:,1)),2:3);
% wellcoords = 
% distanceFromWellAtRewTime = 
% rewTimes = rewTimes(distanceFromWellAtRewTime < 20);

end

function times = get_dios_times(dios, DIO, day, epoch)
o = [];
original_ids = cellfetch(DIO{day}{epoch}, 'original_id');
times = [];
for odio = 1:length(dios)
    o(odio) = find(cell2mat(cellfun(@(x) ...
        strcmp(sprintf('Dout%d', dios(odio)), x), ...
        original_ids.values, 'un', 0)));
    % get rid of artificial first and only keep DIO OUT UP
    t = DIO{day}{epoch}{o(odio)}.times(2:2:end);
    times = [times; t repmat(dios(odio),length(t),1)];
end
if ~isempty(times)
    [~,idx] = sort(times(:,1)); % sort just the first column
    times = times(idx,:);   % sort the whole matrix using the sort indices
end
end