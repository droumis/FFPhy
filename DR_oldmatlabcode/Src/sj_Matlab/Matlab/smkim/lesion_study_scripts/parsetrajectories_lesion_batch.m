load('groups');
subjects = {groups(find([groups(:).keep])).subject};

for s = 1:numel(subjects)
    disp(['processing ' subjects{s}]);
%
    load([subjects{s} '_lineartrack_smoothedpos.mat']);
    for i = 1:2
        for j = 1:2
            try
                journeys{i}{j} = parsetrajectories_lineartrack(smoothedpos{i}{j},0.5,15);
            end
        end
    end
    save([subjects{s} '_lineartrack_journeys.mat'],'journeys');
    clear('journeys','smoothedpos');
%
    load([subjects{s} '_Wtrack_smoothedpos.mat']);
    for i = 1:10
        for j = 1:2
            try
                journeys{i}{j} = parsetrajectories_Wtrack(smoothedpos{i}{j},1,15);
            end
        end
    end
    save([subjects{s} '_Wtrack_journeys.mat'],'journeys');
    clear('journeys','smoothedpos');
%{
    load([subjects{s} '_Wtrack_smoothedpos.mat']);
    for i = 1:10
        for j = 1:2
            try
                passes{i}{j} = parsetrajectories_Wtrack_choicepoint(smoothedpos{i}{j},0.5,5);
            end
        end
    end
    save([subjects{s} '_Wtrack_choicepoint.mat'],'passes');
    clear('passes','smoothedpos');
%}
end
