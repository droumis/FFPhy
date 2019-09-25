%PLOT SPIKEWIDTH VS RATE FOR ALL ANIMALS
for a = 1:11
    switch a
        case 1
            animalinfo = {'Bond','/data13/mcarr/Bon','bon'};
        case 2
            animalinfo = {'Conley','/data13/mcarr/Con','con'};
        case 3
            animalinfo = {'Corriander','/data13/mcarr/Cor','Cor'};
        case 4
            animalinfo = {'Dudley','/data13/mcarr/Dud/','dud'};
        case 5
            animalinfo = {'Eight','/data13/mcarr/Eig/','Eig'};
        case 6
            animalinfo = {'Five','/data13/mcarr/Fiv/','Fiv'};
        case 7
            animalinfo = {'Frank','/data13/mcarr/Fra/','fra'};
        case 8
            animalinfo = {'Miles','/data13/mcarr/Mil/','mil'};
        case 9
            animalinfo = {'Seven','/data13/mcarr/Sev/','Sev'};
        case 10
            animalinfo = {'Six','/data13/mcarr/Six/','Six'};
        case 11
            animalinfo = {'Ten','/data13/mcarr/Ten/','ten'};
    end
    out{a} = getspikewidthvrate(animalinfo{2},animalinfo{3});
end

figure
hold on
for a = 1:length(out)
    plot(out{a}(:,2),out{a}(:,1),'.')
    ylabel('Spike Width');
    xlabel('Mean Rate');
end

% Based on the plot, use separation line to define each cells as fs if they
% exceed mean rate & spike width.
for a = 1:11
    switch a
        case 1
            animalinfo = {'Bond','/data13/mcarr/Bon','bon'};
        case 2
            animalinfo = {'Conley','/data13/mcarr/Con','con'};
        case 3
            animalinfo = {'Corriander','/data13/mcarr/Cor','Cor'};
        case 4
            animalinfo = {'Dudley','/data13/mcarr/Dud/','dud'};
        case 5
            animalinfo = {'Eight','/data13/mcarr/Eig/','Eig'};
        case 6
            animalinfo = {'Five','/data13/mcarr/Fiv/','Fiv'};
        case 7
            animalinfo = {'Frank','/data13/mcarr/Fra/','fra'};
        case 8
            animalinfo = {'Miles','/data13/mcarr/Mil/','mil'};
        case 9
            animalinfo = {'Seven','/data13/mcarr/Sev/','Sev'};
        case 10
            animalinfo = {'Six','/data13/mcarr/Six/','Six'};
        case 11
            animalinfo = {'Ten','/data13/mcarr/Ten/','ten'};
    end

load([animalinfo{2}, '/', animalinfo{3},'cellinfo']);

meanrate = cellfetch(cellinfo,'meanrate');
spikewidth = cellfetch(cellinfo,'spikewidth');
index = meanrate.index;
for j = 1:length(index)
    if ~isempty(meanrate.values{j}) && ~isempty(spikewidth.values{j})
        fs = meanrate.values{j}>5 & spikewidth.values{j} < 10;
        if fs
            cellinfo{index(j,1)}{index(j,2)}{index(j,3)}{index(j,4)}.fs = 1;
        end
    end
end

save([animalinfo{2},animalinfo{3},'cellinfo'],'cellinfo');
end