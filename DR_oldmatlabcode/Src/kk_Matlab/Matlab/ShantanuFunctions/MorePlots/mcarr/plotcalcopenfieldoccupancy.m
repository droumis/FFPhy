%Animal selection
%-----------------------------------------------------
animals = {'Eight'};
%-----------------------------------------------------


%Filter creation%
%--------------------------------------------------------
epochfilter = '(isequal($description, ''TrackA''))';

cellfilter = '(isequal($area, ''MEC''))';

timefilter = { {'getlinstate', '(($traj ~= -1) & abs($velocity) >= 1)', 6} };

cell1f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

%-----------------------------------------------------------

%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

cell1f = setfilteriterator(cell1f,iterator);

%non occ-norm mean rate
cell1f = setfilterfunction(cell1f, 'calcopenfieldoccupancy', {'spikes', 'pos'},'std', 1.5, 'appendindex', 1);

cell1f = runfilter(cell1f);

%-------------------------------------------------

%%
f = cell1f;
n = 1;

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
        cells = unique(allcells(:,[1 2 3 4]), 'rows');
        for c = 1:size(cells,1)
            cell = cells(c,:); %d e t c
                figure(c)
                subplot(length(f(an).output),1,g)
                imagesc(f(an).output{g}(c).smoothedspikerate, [-0.1 max(max(f(an).output{g}(c).smoothedspikerate))*.75] )
                colormap(cmap)
                xlim([-10 250])
                ylim([-10 100])
                colorbar
                figure(c)
                title([f(an).animal{3}, ' Day', num2str(cell(1)), ' Epochs', num2str(cell(2)), ' Tetrode', num2str(cell(3)), ' Cell', num2str(cell(4))])
            %pause
        end
    end
end