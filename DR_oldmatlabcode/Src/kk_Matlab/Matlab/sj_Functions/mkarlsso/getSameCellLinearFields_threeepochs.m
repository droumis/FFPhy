function out = getSameCellLinearFields_threeepochs(datafilter)
% out = getSameCellLinearFields(datafilter)
% Assumes that datafilter(animal).out has two cells.  Inside each cell is a
% structure of length NUMCELLS with the fields index ([day epoch tetrode cell]) and trajdata.
% This structure is the output of FILTERCALCLINFIELDS.  It is assumed that
% the data in the first cell is for the first run and the data in the
% second cell is for the second run in the pair, and that these runs
% happend on the same day so that two runs can be paired for each cell. OUT
% is a structure with fields RUNINDEX (2 by 4, contains the indices for the
% two runs), trajnum, dist, run1rates, and run2rates.


paircount = 0;
out = struct;
for i = 1:length(datafilter)
    
    if (length(datafilter(i).output) > 1)
        
        %compile both index lists
        g1index = [];
        g2index = [];
        g3index = [];
        for group1index = 1:length(datafilter(i).output{1})
            g1index = stack(g1index,[group1index datafilter(i).output{1}(group1index).index]);
        end
        g1index = g1index(find(~isnan(sum(g1index,2))),:);
        
        for group2index = 1:length(datafilter(i).output{2})
            g2index = stack(g2index,[group2index datafilter(i).output{2}(group2index).index]);
        end
        g2index = g2index(find(~isnan(sum(g2index,2))),:);
        
        for group3index = 1:length(datafilter(i).output{3})
            g3index = stack(g3index,[group3index datafilter(i).output{3}(group3index).index]);
        end
        g3index = g3index(find(~isnan(sum(g3index,2))),:);
        
        
        %find all matching index pairs (ignores the epoch)
        matches = rowfind(g1index(:,[2 4 5]),g2index(:,[2 4 5]));
        matches2 = rowfind(g1index(:,[2 4 5]),g3index(:,[2 4 5]));
        
        
        %go through each match and compile the data from the two runs
        for j = 1:length(matches)
            if ((matches(j) > 0) & (matches2(j) > 0))
                paircount = paircount+1;
                trajnum = [];
                dist = [];
                run1rates = [];
                run2rates = [];
                run3rates = [];
                %get the indices for the two runs
                runindex = [ [i g1index(j,2:5)]; [i g2index(matches(j),2:5)]; [i g3index(matches2(j),2:5)]];
                trajdata1 = datafilter(i).output{1}(g1index(j,1)).trajdata; %get the traj data in the datafilter for pair
                trajdata2 = datafilter(i).output{2}(g2index(matches(j),1)).trajdata;
                trajdata3 = datafilter(i).output{3}(g3index(matches2(j),1)).trajdata;
                for traj = 1:length(trajdata1)
                    
                    %if the lengths are different, we crop the shorter one
                    trajlengths = [length(trajdata1{traj}(:,5)) length(trajdata2{traj}(:,5)) length(trajdata3{traj}(:,5))];

                    trajlength = min(trajlengths);
                    tmptrajnum = ones(trajlength,1)*traj;
                    trajnum = [trajnum;tmptrajnum];

                    dist = [dist; trajdata1{traj}(1:trajlength,1)];

                    run1rates = [run1rates; trajdata1{traj}(1:trajlength,5)];
                    run2rates = [run2rates; trajdata2{traj}(1:trajlength,5)];
                    run3rates = [run3rates; trajdata3{traj}(1:trajlength,5)];

                end
                out(paircount).runindex = runindex;
                out(paircount).trajnum = trajnum;
                out(paircount).dist = dist;
                out(paircount).run1rates = run1rates;
                out(paircount).run2rates = run2rates;
                out(paircount).run3rates = run3rates;
            end
        end
    end
    
end