function rewTimes = get_rewardTimes(task, DIO, day, epochs)


outdios = task{day}{epochs(1)}.outputdio;
original_ids = cellfetch(DIO{day}{epochs(1)}, 'original_id');
rewTimes = [];
o = [];
for e = 1:length(epochs)
    for odio = 1:length(outdios)
        o(odio) = find(cell2mat(cellfun(@(x) ...
            strcmp(sprintf('Dout%d', outdios(odio)),x), ...
            original_ids.values, 'un', 0)));
        rewTimes = [rewTimes; DIO{day}{epochs(e)}{o(odio)}.times ...
            repmat(outdios(odio),length(DIO{day}{epochs(e)}{o(odio)}.times)...
            ,1)];
    end
end
% get rid of duplicate timestamps in output.. (trodes bug)
[R, uid, ~] = unique(rewTimes(:,1));
R = [R rewTimes(uid,2)];
[~, sidx] = sort(R(:,1));
rewTimes = R(sidx,:);
rewTimes = rewTimes(2:end,:); % drop artificial first

% confirm that the animal is at the reward well for each rewTime
% condition: animal pos is within X of well at time of rew
% condition: there is a well input trigger 1 second preceding rew time
% rewTimes = rewTimes(distanceFromWellAtRewTime < 20);
% posAtRewTime = pos(lookup(rewTimes(:,1), pos(:,1)),2:3);
% wellcoords = 
% distanceFromWellAtRewTime = 
% rewTimes = rewTimes(distanceFromWellAtRewTime < 20);

end