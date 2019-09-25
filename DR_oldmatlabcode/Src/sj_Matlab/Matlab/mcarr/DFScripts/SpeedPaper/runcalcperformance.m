%Animal Selection
animals = {'Conley','Corriander','Dudley','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];

for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Time Filter
timefilter = [];

% Behavior Filter

iterator = 'singleepochanal';

f = createfilter('animal',animals,'epochs',epochfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcperformance', {'linpos','task'}, [2 1 3],'includetrajbound',1);

f = runfilter(f);

% PLOT
tmpcorrect = numericgroupcombine(f);
correct = zeros(size(tmpcorrect));
se_correct = zeros(size(correct));
for i = 1:length(correct)
    correct(i) = mean(1 - tmpcorrect{i}(:,2));
    se_correct(i) = std(1 -tmpcorrect{i}(:,2))./sqrt(size(tmpcorrect{i},1) + 1);
end
errorU = [se_correct + correct];
errorL = [correct - se_correct];
error = [errorU fliplr(errorL)];

figure
plot(1:length(correct),correct,'k')
hold on
fill([1:length(correct) fliplr(1:length(correct))], error,'k')
box off
set(gca,'xlim',[0 11],'ylim', [0 0.5],'xtick',[1:10])
xlabel('Exposure')
ylabel('Proportion Error Trials')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/behavioralperformance.pdf', m, d, y);
print('-dpdf', savestring)
