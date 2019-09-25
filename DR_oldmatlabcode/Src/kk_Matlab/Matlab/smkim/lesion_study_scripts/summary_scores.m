

lesion_subjects = { 'M02' 'M03' 'M12' 'M13' 'M14' 'M16' 'M17' 'M19' 'M20' 'M22' };
control_subjects= { 'M06' 'M24' 'M25' 'M26' };

% 10 days, 2 sessions per day

lesion_correct = NaN(10,numel(lesion_subjects));
control_correct = NaN(10,numel(control_subjects));

for i = 1:length(lesion_subjects)
    currentsubject = lesion_subjects{i};
    currentscores = [ currentsubject '_Wtrack_scores.mat'];
    load(currentscores);
    for j = 1:10
        session1_scores = scores{j}{1}.data(find( (scores{j}{1}.data(:,1)==2) & (scores{j}{1}.data(:,2)~=2) ),3);
        session2_scores = scores{j}{2}.data(find( (scores{j}{2}.data(:,1)==2) & (scores{j}{2}.data(:,2)~=2) ),3);

%        session1_scores = scores{j}{1}.data(find(scores{j}{1}.data(:,1)~=2),3);
%        session2_scores = scores{j}{2}.data(find(scores{j}{2}.data(:,1)~=2),3);

%        session1_scores = scores{j}{1}.data(:,3);
%        session2_scores = scores{j}{2}.data(:,3);

        day_scores = [session1_scores; session2_scores];
        lesion_correct(j,i) = nnz(day_scores)/numel(day_scores);
    end
    clear(currentscores);
end

for i = 1:length(control_subjects)
    currentsubject = control_subjects{i};
    currentscores = [ currentsubject '_Wtrack_scores.mat'];
    load(currentscores);
    for j = 1:10
        session1_scores = scores{j}{1}.data(find( (scores{j}{1}.data(:,1)==2) & (scores{j}{1}.data(:,2)~=2) ),3);
        session2_scores = scores{j}{2}.data(find( (scores{j}{2}.data(:,1)==2) & (scores{j}{2}.data(:,2)~=2) ),3);

%        session1_scores = scores{j}{1}.data(find(scores{j}{1}.data(:,1)~=2),3);
%        session2_scores = scores{j}{2}.data(find(scores{j}{2}.data(:,1)~=2),3);

%        session1_scores = scores{j}{1}.data(:,3);
%        session2_scores = scores{j}{2}.data(:,3);

        day_scores = [session1_scores; session2_scores];
        control_correct(j,i) = nnz(day_scores)/numel(day_scores);
    end
end

lesion_mean = mean(lesion_correct,2);
lesion_stdev = std(lesion_correct,0,2);
lesion_sem = lesion_stdev/sqrt(size(lesion_correct,2));

control_mean = mean(control_correct,2);
control_stdev = std(control_correct,0,2);
control_sem = control_stdev/sqrt(size(control_correct,2));

% plot summary
figure(1);
set(gcf,'Color','w');
set(gca,'FontSize',20,'LineWidth',3,'Color','w','XLim',[0.5 10.5],'YLim',[0 1], ...
    'XColor','k','YColor','k','XTick',1:10);
hold on;

control_lines = errorbar(1:10,control_mean,control_sem, ...
    'Color','k','LineStyle','-','Marker','s','MarkerFaceColor','k','LineWidth',3,'MarkerSize',12);
lesion_lines = errorbar(1:10,lesion_mean,lesion_sem, ...
    'Color','k','LineStyle','--','Marker','o','MarkerFaceColor','w','LineWidth',3,'MarkerSize',12);

set(get(gca,'XLabel'),'FontSize',20,'Color','k','String','day');
set(get(gca,'YLabel'),'FontSize',20,'Color','k','String','proportion correct food well visits');

legend_handle = legend('sham-operated control (n=4)','lesion (n=10)', ...
    'Location','SouthEast');
set(legend_handle,'Box','off','TextColor','k');

