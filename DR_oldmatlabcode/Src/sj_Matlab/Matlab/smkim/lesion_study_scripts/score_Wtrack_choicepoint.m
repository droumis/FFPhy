load('Wtrack_summary');

for i = 1:length(Wtrack_summary)
    subject = Wtrack_summary(i).subject;
    load([subject '_Wtrack_choicepoint.mat']);
    load([subject '_Wtrack_trajectories.mat']);
    for j = 1:10
        t.startzone = [trajectories{j}{1}.startzone, trajectories{j}{2}.startzone];
        t.endzone = [trajectories{j}{1}.endzone, trajectories{j}{2}.endzone];
        Wtrack_summary(i).revisit(j) = nnz(t.startzone == t.endzone) / Wtrack_summary(i).total(j);
        Wtrack_summary(i).right(j) = nnz(t.startzone == 1);
        Wtrack_summary(i).center(j) = nnz(t.startzone == 2);
        Wtrack_summary(i).left(j) = nnz(t.startzone == 3);
        p.speed = [ passes{j}{1}.speed, passes{j}{2}.speed ];
        p.from_arm = [ passes{j}{1}.from_arm; passes{j}{2}.from_arm ];
        p.to_arm = [ passes{j}{1}.to_arm; passes{j}{2}.to_arm ];
        choicepoint.total(j) = length(p.from_arm);
        choicepoint.fromside(j) = nnz(p.from_arm ~= 2);
        choicepoint.fromside_turnback(j) = nnz( ...
            (p.from_arm ~= 2) & (p.from_arm == p.to_arm) );
        choicepoint.fromside_skipcenter(j) = nnz( ...
            (p.from_arm ~= 2) & (p.to_arm ~= 2) & (p.from_arm ~= p.to_arm) );
        choicepoint.speed(j) = mean(p.speed);
    end
    Wtrack_summary(i).choicepoint = choicepoint;
    clear('passes');
end
%save([subject '_Wtrack_summary.mat'],'Wtrack_summary');

