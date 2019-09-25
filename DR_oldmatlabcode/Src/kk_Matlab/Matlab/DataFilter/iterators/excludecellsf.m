function f = excludecellsf(f, includecells)
% only includes cells in includecells, can apply to single cells or cell
% pairs

%iterate through all animals
for an = 1:length(f)

    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)

        for e = 1:size(f(an).epochs{g},1) %for each epoch
            day = f(an).epochs{g}(e,1);
            epoch = f(an).epochs{g}(e,2);
            includetetcell = includecells(includecells(:,1)==an & includecells(:,2)==day & includecells(:,3)==epoch, 4:5);
            tmpdata = f(an).data{g}{e};
            if size(tmpdata, 2) ==2
                newtmpdata = tmpdata(ismember(tmpdata, includetetcell, 'rows'),:);
            elseif size(tmpdata, 2) ==4
                newtmpdata = tmpdata(ismember(tmpdata(:,1:2), includetetcell, 'rows') & ismember(tmpdata(:,3:4), includetetcell, 'rows'), :);
            end
           f(an).data{g}{e} = newtmpdata;
        end
    end
end



