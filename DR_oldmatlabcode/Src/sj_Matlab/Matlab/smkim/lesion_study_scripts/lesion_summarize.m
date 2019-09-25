load('groups');

subjects = {groups(find([groups(:).keep])).subject};

for s = 1:length(subjects)
    %lineartrack_summary: subject, group, total, correct, left, right, speed
    load([subjects{s} '_lineartrack_journeys.mat']);
    lineartrack_summary(s).subject = subjects{s};
    lineartrack_summary(s).group = groups(find(strcmp({groups(:).subject},subjects{s}))).group;
    for d = 1:2
        lineartrack_summary(s).total(d) = ...
            nnz(~isnan(journeys{d}{1}.correct)) + nnz(~isnan(journeys{d}{2}.correct));
        lineartrack_summary(s).correct(d) = ...
            nnz(journeys{d}{1}.correct == 1) + nnz(journeys{d}{2}.correct == 1);
        lineartrack_summary(s).left(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & (journeys{d}{1}.endzone == 2)) + ... 
            nnz(~isnan(journeys{d}{2}.correct) & (journeys{d}{2}.endzone == 2));
        lineartrack_summary(s).right(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & (journeys{d}{1}.endzone == 1)) + ... 
            nnz(~isnan(journeys{d}{2}.correct) & (journeys{d}{2}.endzone == 1));
        lineartrack_summary(s).perjourneyspeed(d) = median([ ...
            journeys{d}{1}.speed(~isnan(journeys{d}{1}.correct)), ...
            journeys{d}{2}.speed(~isnan(journeys{d}{2}.correct)) ]);
        concat_snippets = [ ...
            vertcat(journeys{d}{1}.smoothedpos_snippets(:).data); ...
            vertcat(journeys{d}{2}.smoothedpos_snippets(:).data) ];
        lineartrack_summary(s).overallspeed(d) = mean(hypot( ...
            concat_snippets(:,4), concat_snippets(:,5) ));
        totalruntime = ...
            timetrans({journeys{d}{1}.tend},1,2)-timetrans({journeys{d}{1}.tstart},1,2) + ...
            timetrans({journeys{d}{2}.tend},1,2)-timetrans({journeys{d}{2}.tstart},1,2);
        timeatfoodwells = totalruntime - size(concat_snippets,1)/30;
        lineartrack_summary(s).foodwelloccupancy(d) = timeatfoodwells/totalruntime;
        lineartrack_summary(s).visitduration(d) = median([ ...
            journeys{d}{1}.visitduration(~isnan(journeys{d}{1}.visitduration)), ...
            journeys{d}{2}.visitduration(~isnan(journeys{d}{2}.visitduration)) ]);
    end
        
    %Wtrack_summary: subject, group, total, correct, left, right, center, revisit, choicepoint, speed
    load([subjects{s} '_Wtrack_journeys.mat']);
    Wtrack_summary(s).subject = subjects{s};
    Wtrack_summary(s).group = groups(find(strcmp({groups(:).subject},subjects{s}))).group;
    for d = 1:10
        Wtrack_summary(s).total(d) = ...
            nnz(~isnan(journeys{d}{1}.correct)) + nnz(~isnan(journeys{d}{2}.correct));
        Wtrack_summary(s).correct(d) = ...
            nnz(journeys{d}{1}.correct == 1) + nnz(journeys{d}{2}.correct == 1);
        Wtrack_summary(s).left(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & (journeys{d}{1}.endzone == 3)) + ... 
            nnz(~isnan(journeys{d}{2}.correct) & (journeys{d}{2}.endzone == 3));
        Wtrack_summary(s).right(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & (journeys{d}{1}.endzone == 1)) + ... 
            nnz(~isnan(journeys{d}{2}.correct) & (journeys{d}{2}.endzone == 1));
        Wtrack_summary(s).center(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & (journeys{d}{1}.endzone == 2)) + ... 
            nnz(~isnan(journeys{d}{2}.correct) & (journeys{d}{2}.endzone == 2));
        Wtrack_summary(s).revisit(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & ...
            (journeys{d}{1}.endzone == journeys{d}{1}.startzone)) + ... 
            nnz(~isnan(journeys{d}{2}.correct) & ...
            (journeys{d}{2}.endzone == journeys{d}{2}.startzone));
        Wtrack_summary(s).choicepoint(d) = ...
            nnz(~isnan(journeys{d}{1}.correct) & journeys{d}{1}.choicepoint) + ...
            nnz(~isnan(journeys{d}{2}.correct) & journeys{d}{2}.choicepoint);
        Wtrack_summary(s).perjourneyspeed(d) = median([ ...
            journeys{d}{1}.speed(~isnan(journeys{d}{1}.correct)), ...
            journeys{d}{2}.speed(~isnan(journeys{d}{2}.correct)) ]);
        concat_snippets = [ ...
            vertcat(journeys{d}{1}.smoothedpos_snippets(:).data); ...
            vertcat(journeys{d}{2}.smoothedpos_snippets(:).data) ];
        Wtrack_summary(s).overallspeed(d) = mean(hypot( ...
            concat_snippets(:,4), concat_snippets(:,5) ));
        totalruntime = ...
            timetrans({journeys{d}{1}.tend},1,2)-timetrans({journeys{d}{1}.tstart},1,2) + ...
            timetrans({journeys{d}{2}.tend},1,2)-timetrans({journeys{d}{2}.tstart},1,2);
        timeatfoodwells = totalruntime - size(concat_snippets,1)/30;
        Wtrack_summary(s).foodwelloccupancy(d) = timeatfoodwells/totalruntime;
        Wtrack_summary(s).visitduration(d) = median([ ...
            journeys{d}{1}.visitduration(~isnan(journeys{d}{1}.visitduration)), ...
            journeys{d}{2}.visitduration(~isnan(journeys{d}{2}.visitduration)) ]);
    end
end    
save('lineartrack_summary','lineartrack_summary');
save('Wtrack_summary','Wtrack_summary');
