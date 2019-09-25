Veqn = '>=3'
minVPF = 2 %cm/sec
minPeakPF = 3
mintime = 10
n = 1


%Animal selection
%-----------------------------------------------------
%animals = {'Dwight'} %'Barack', 'Calvin', 'Dwight'};
%animals = {'Miles'}%, 'Ten', 'Conley', 'Dudley'}; %Mattias's animal: 1st exposure to first track on first day
%animals = {'Miles'};
animals = {'Dudley'};
%-----------------------------------------------------


%Filter creation%
%--------------------------------------------------------
epochfilter{1} = ['(isequal($type, ''run''))']; %just analyze days where switching between tasks

cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 4))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal

timefilter = { {'getlinvelocity', ['((abs($velocity)', Veqn,'))']} };

cell1f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
cell1f = testexcludetimes(cell1f, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30

%only include cells with placefields
if minPeakPF>0
    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfilter);
    cell1f = excludecellsf(cell1f, includecells);
end
%-----------------------------------------------------------

%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

cell1f = setfilteriterator(cell1f,iterator);

%non occ-norm mean rate
cell1f = setfilterfunction(cell1f, 'calcopenfieldoccupancy', {'spikes', 'pos'},'std', 1.5, 'appendindex', 1);

cell1f = runfilter(cell1f);


%cellgroups = numericgroupcombine(cell1f, 1);  %if output numeric, this combines across animals, first col is animal
%-------------------------------------------------

%%
f = cell1f;

colormap('default');
cmap = jet(1024) ./ 1.5;
cmap = cmap(100:920,:); 
cmap(1,:) = 1;
colormap(cmap);

for an = 1:length(f) %for each animal
    for g = 1:length(f(an).output) %for each epoch grouping/condition
        %find all identical cells in diff epochs
        allcells = [];
        for e = 1:size(f(an).epochs{g},1) %for each epoch
            de = f(an).epochs{g}(e,:); %day epoch
            allcells = [allcells; repmat(de, size(f(an).data{g}{e},1),1) f(an).data{g}{e}]; %d e t c for all cells in g, indexes should correspond to indexes of f(an).output{g}
        end
        cells = unique(allcells(:,[1 3 4]), 'rows'); %remove repeats due to epochs
        for c = 1:size(cells,1)
            cell = cells(c,:); %d t c
            ind = find(ismember(allcells(:,[1 3 4]), cell,'rows'));
            for i = 1:length(ind)
                figure(n)
                subplot(length(ind), 1,i)
                imagesc(f(an).output{g}(ind(i)).smoothedspikerate, [-0.1 max(max(f(an).output{g}(ind(i)).smoothedspikerate))*.75] )
                colormap(cmap)
                xlim([-10 250])
                ylim([-10 100])
                colorbar
                figure(  n+1)
                subplot(length(ind), 1,i)
                imagesc(f(an).output{g}(ind(i)).smoothedspikerate, [-0.1 3])
                colormap(cmap)
                xlim([-10 250])
                ylim([-10 100])
                colorbar
            end
            epochs = allcells(ind,2);
            figure(n)
            subtitle([f(an).animal{3}, ' Day', num2str(cell(1)), ' Epochs', num2str(epochs'), ' Tetrode', num2str(cell(2)), ' Cell', num2str(cell(3))])
            figure(n+1)
            subtitle([f(an).animal{3},' Day', num2str(cell(1)), ' Epochs', num2str(epochs'), ' Tetrode', num2str(cell(2)), ' Cell', num2str(cell(3))])
            pause

            figure(n)
            clf
            figure(n+1)
            clf
        end
    end
end

%% print to ps
%deinfe new cells
if exist('newcells')
for c = 1:size(newcells,1)
    cell = newcells(c,:); %d t c
    ind = find(ismember(allcells(:,[1 3 4]), cell,'rows'));
    for i = 1:length(ind)
        figure(n)
        subplot(length(ind), 1,i)
        imagesc(f(an).output{g}(ind(i)).smoothedspikerate, [-0.1 max(max(f(an).output{g}(ind(i)).smoothedspikerate))*.75] )
        colormap(cmap)
        xlim([-10 250])
        ylim([-10 100])
        colorbar
        
        figure(n+1)
        subplot(length(ind), 1,i)
        imagesc(f(an).output{g}(ind(i)).smoothedspikerate, [-0.1 3])
        colormap(cmap)
        xlim([-10 250])
        ylim([-10 100])
        colorbar
    end
    epochs = allcells(ind,2);
    figure(n)
    titl = [f(an).animal{3}, 'D', num2str(cell(1)), 'T', num2str(cell(2)), 'C', num2str(cell(3))];
    subtitle(titl)
    print('-dpsc', titl)
    figure(n+1)
    titl = [f(an).animal{3}, 'D', num2str(cell(1)), 'T', num2str(cell(2)), 'C', num2str(cell(3)), '_3HzThresh'];
    subtitle(titl)
    print('-dpsc', titl)
    n=n+2;
end
end

