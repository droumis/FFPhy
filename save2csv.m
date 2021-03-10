



animal = 'JZ2';
andef = animaldef(animal);

%% load behave state and convert mat to table
bstate = load_data('filterframework', [], animal, 'filtfunction', 'behavestate');
T = array2table(bstate.statespace.allepsMat);
H = bstate.statespace.allepsMatFields;

%% add var headers
T.Properties.VariableNames = strsplit(H,' ');
T.animal(:,1) = repmat(string(animal), height(T), 1);

%% create day and epoch vecs
dayeps = bstate.statespace.allbound(:,[5 6]);
% eplens = bstate.statespace.eplengths;
% days = 1:(length(eplens) / 2);
% eps = repmat([2;4], (length(eplens) / 2), 1);
%%
taskInfo = {};
days = sort(unique(dayeps(:,1)));
for d = 1:length(days)
    itask = load_data('filterframework', sprintf('task%02d',days(d)), animal, 'animpos', 0);
    taskInfo{d} = itask.task{d};
end
%%
% dayscat = [days; days];
% dayscat = dayscat(:);
% 
% alldays = [];
% alleps  = [];
% allenvs = [];
% alldaytype = [];
% 
% for e = 1:length(bstate.statespace.eplengths)
%     elen = bstate.statespace.eplengths(e);
%     alldays = [alldays ; repmat(dayscat(e), elen, 1)];
%     alleps = [alleps ; repmat(eps(e), elen, 1)];
% %     allenvs = [allenvs; repmat(string(taskInfo{dayscat(e)}{eps(e)}.environment), elen, 1)];
% %     alldaytype = [alldaytype; repmat(taskInfo{dayscat(e)}{eps(e)}.daytype, elen, 1)];
% end
%%
T.day = dayeps(:,1);
T.epoch = dayeps(:,2);
%%
for i = 1:height(T)
%    fprintf('%d\n',i)
    try
        T.environment(i) = string(taskInfo{T.day(i)}{T.epoch(i)}.daytype);
    catch
        if strcmp(animal, 'D12') && T.day(i) == 7
            T.environment(i) = 'wtrack';
        end
    end
end
%%

writetable(T,sprintf('/home/droumis/Src/ephys_autoW/Processed-Data/%s/csv/%s_wellVisits.csv', animal, animal))
