Veqn = '>3'
minV=str2num(Veqn(end))
minVPF = 3 %cm/sec
minPeakPF = 3
mintime = 3
n = 1
inclStates = 6
traj = [3 4]
wtrack = 0; % 
%traj = [1 4; 2 3] %[1:4]; %[1 4; 2 3];
%armlength = 66 %70 Wtrack %66 6 arm maze;
maxstage = [12 3]
comp = 1
correct = 0
excluderepeat = 1
area = 'CA3'
fnc = 3

%Animal selection
%-----------------------------------------------------
animals = {'Barack', 'Calvin', 'Dwight'};
if wtrack ==1
    animals = {'Bond', 'Frank', 'Nine'} %for more familiar track A
    %animals = {'Miles', 'Ten', 'Conley', 'Dudley'}; %Mattias's animal: 1st exposure to first track on first day
end
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
for i = maxstage;
    epochfilter{1} = ['($switchday >= 1)']; %just analyze days where switching between tasks
    %epochfilter{1} = ['($exposure >= 1)  & ($exposure <= 3)  & ($tasknum == 1 )']; %just analyze days where switching between tasks
    epochfPF = ['($switchday > -6)']; %just analyze days where switching between tasks
    
    if wtrack ==1
        epochfilter{1} = ['(isequal($type, ''run''))  & (isequal($environment, ''TrackB''))']; %just analyze days where switching between tasks%
        epochfPF = epochfilter{1} ;
    end
    
    if isequal(area, 'CA1')
        cellfilter = '(isequal($area, ''CA1'') & ($meanrate < 4))'
    else
        cellfilter = '(isequal($area, ''CA3'') & ($meanrate < 4))'
    end
    
    timefilter =  {{'getdistclosestwell', '($distwell >= 10)'} {'getlinvelocity', ['((abs($velocity) ',Veqn,'))'] } };
    if correct == 1 && wtrack == 0
        timefilter{length(timefilter)+1} = {'filtergetcorrecttraj', '($state > 0)', 6, minV, 0, 10, 150} ;
    elseif correct == 1 && wtrack == 1
        timefilter{length(timefilter)+1} = {'filtergetcorrecttraj', '($state > 0)', 6, minV, 0, 10, 150, 'correctorder', [2 1 3]} ;
    end
    if i < 10
        timefilter{length(timefilter)+1} = {'getcalctaskstage', ['($includebehave == ', num2str(i), ')'], comp} ;
    elseif i == 12
        timefilter{length(timefilter)+1} = {'getcalctaskstage', ['($includebehave == 1) | ($includebehave == 2)'], comp} ;
    end
    
    
    cell1f{i} = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
    cell1f{i} = testexcludetimes(cell1f{i}, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30
    
    %only include cells with placefields
    if minPeakPF>0
        lessthan = 0;
        includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan, area);
        if excluderepeat == 1
            [b row col] = unique(includecells(:,[1 2 4 5]),'rows', 'first');
            includecells = includecells(row,:);
        end
        
        cell1f{i} = excludecellsf(cell1f{i}, includecells);
    end
    %-----------------------------------------------------------
    
    
    %Run function- single cells
    %--------------------------------------------
    iterator = 'singlecellanal';
    cell1f{i} = setfilteriterator(cell1f{i},iterator);
    
    
    
    
    
    switch (fnc)
            case 1
                %non occ-norm mean rate
                cell1f{i} = setfilterfunction(cell1f{i}, 'calctotalmeanrate', {'spikes'});

            case 2
                %occ-norm mean rate
                cell1f{i} = setfilterfunction(cell1f{i}, 'calcoccnormmeanrate', {'spikes', 'linpos'});

            case 3
                % %peak rate
                cell1f{i} = setfilterfunction(cell1f{i}, 'calcpeakrate2', {'spikes', 'linpos'}); %already have a calcpeakrate function that is not compatible with this anal --> renamed mattias' function calcpeakrate2

        case 4
            cell1f{i} = setfilterfunction(cell1f{i}, 'getspatialinfo', {'linpos', 'spikes'},6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', traj);
        end
    
     cell1f{i} = runfilter(cell1f{i});
end

%% group
for i = maxstage;
    cellgroups{i} = numericgroupcombine(cell1f{i}, 1);  %if output numeric, this combines across animals, first col is animal
end
%% plot

if fnc == 4
data1 = cellgroups{12}{1}(ismember(cellgroups{12}{1}(:,6), 1:6),end);
data2 = cellgroups{3}{1}(ismember(cellgroups{3}{1}(:,6), 1:6),end);

else
data1 = cellgroups{12}{1}(:,end);
data2 = cellgroups{3}{1}(:,end);
end
edges = linspace(min([data1; data2]),max([data1; data2]),10);
h1 = histc(data1,edges);
h2 = histc(data2,edges);
figure
hold on
plot(edges, h1/sum(h1), 'b', 'linewidth', 3)
plot(edges, h2/sum(h2), 'r', 'linewidth', 3)
xlabel(cell1f{12}(1).function.name)
ylabel('fraction')
legend('stage 1 & 2', 'stage 3')

[h p] = kstest2(data1(~isnan(data1)),data2(~isnan(data2)))
med12 = median(data1(~isnan(data1)))
med3 = median(data2(~isnan(data2)))
sizes = [size(data1(~isnan(data1)),1) size(data2(~isnan(data2)),1) ]

rp = plotmeanse(data1, data2, 'stg1&2', 'stg3', cell1f{12}(1).function.name)


