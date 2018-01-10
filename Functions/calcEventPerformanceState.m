
function [out] = calcEventPerformanceState(F, animals, behavestruct)

out.eventstate = [];
out.animals = [];
out.index = [];
disp('CALCULATING EVENTSTATE ++++++++++');
for ianimal = 1:length(F)
    animalinfo = animaldef(lower(animals{ianimal}));
    animalID = animalinfo{1,3}; %use anim prefix for name
    out.animals = [out.animals; {animalID}];
    load(sprintf('%s%s%s.mat',animalinfo{1,2}, animalID, behavestruct));
    
    eventstate = [];
    eventstate.eventStartTimes = [];
    eventstate.state = [];
    eventstate.fields = [];
    out.index{ianimal} = [];
    %get the animal state for every event
    andays = find(~cellfun(@isempty,F(ianimal).output)); %get nonempty eps
    for iday = 1:length(andays)
        day = andays(iday);
        eps = find(~cellfun(@isempty,{F(ianimal).output{day}.index})); %get nonempty eps
        % cross check for wtrack epochs
        load(sprintf('%s%s%s%02d.mat',animalinfo{1,2}, animalID, 'task', day));
        taskinfoall = cellfetch(task, 'environment', 'alltags', 0);
        idayeps = [repmat(day, length(eps),1) eps'];
        tagIndMap = [];
        [~, tagIndMap] = ismember(idayeps,taskinfoall.index(:,[1:2]), 'rows');
        if ~any(tagIndMap)
            continue
        end
        taskTags = taskinfoall.values(tagIndMap);
        useeps = find(cell2mat(strfind(taskTags, 'wtrack')));
        if isempty(useeps)
            continue
        else
            eps = eps(useeps);
        end
        out.index{ianimal} = [out.index{ianimal}; F(ianimal).output{day}(eps(1)).index];
        %load linpos per day
        load(sprintf('%s%s%s%02d.mat',animalinfo{1,2}, animalID, 'linpos', day));
        for iep = eps;
            % get trajectory, segment, and linddist
            eventStartIndices = F(ianimal).output{day}(iep).eventStartIndices;
            LFPtimes = F(ianimal).output{day}(iep).LFPtimes;
            eventStartTimes = LFPtimes(eventStartIndices);
            statematrix = linpos{day}{iep}.statematrix;
            statemat.data = [statematrix.time statematrix.traj statematrix.segmentIndex statematrix.lindist];
            statemat.fields = 'postime traj segment lindist';
            stateindex = knnsearch(statemat.data(:,1), eventStartTimes);
            eventTrajs = statemat.data(stateindex,:);
            %return corresponding vectors for 'correct' and 'mistake' performance state for given day, ep
            [corrvec, mistvec] = getPerformanceState(BehaveState, day, iep, LFPtimes);
            eventindex = knnsearch(LFPtimes, eventStartTimes);
            tmpeventstate = [];
            tmpeventstate = [LFPtimes(eventindex) eventTrajs corrvec(eventindex) mistvec(eventindex) repmat([day iep], length(eventindex),1)];
            eventstate.state = [eventstate.state; tmpeventstate];
            eventstate.fields = {'LFPtime', 'postime', 'traj', 'segment', 'lindist', 'correct', 'mistake', 'day', 'epoch'};
            eventstate.eventStartTimes = [eventstate.eventStartTimes; eventStartTimes];
        end
    end
    if isempty(eventstate.state)
        continue
    end
%     removevec = ones(length(eventstate.state(:,1)),1);
%     if max(abs(diff(eventstate.state(:,1:2),[],2))) > 0.033;
%         disp(sprintf('%d max event-time offset between pos and lfp times is more than 33ms (1 cam frame).. removing offset events', max(abs(diff(eventstate.state(:,1:2),[],2)))))
%         removeevents = find(abs(diff(eventstate.state(:,1:2),[],2)) > 0.033);
% %         removevec = ones(length(ixpc.eventstate{ian}.state(:,1)),1);
%         removevec(removeevents) = 0;
%         eventstate.state = eventstate.state(logical(removevec),:);
%         eventstate.eventStartTimes = eventstate.eventStartTimes(logical(removevec),:);
%     end

%     out.removevec{ianimal} = removevec;
    out.eventstate{ianimal} = eventstate;

end