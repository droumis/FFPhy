function [out] = gettetmaxcell(cellinfo, task, day, cellfilter, includeepochs)

%   includeepochs, 1, 2, 3, default 1
%       1 -- all epochs
%       2 -- all run
%       3 -- all sleeps

% called by getriptimes and getripactivprob and gettetrodelist if using
% maxcell option

%find tetrode with max number of cells per day
if includeepochs == 1
    numcells = [0 0];
    d = day;
    for e = 1:length(cellinfo{d})
        tetlist =  evaluatefilter(cellinfo{d}{e}, cellfilter);
        % get rid of the cell indeces and extract only the tetrode numbers
        if ~isempty(tetlist)
            tetlist = unique(tetlist(:,1))';
            for s = 1:length(tetlist)
                tet = tetlist(s); %tetrode

                if ismember(tet,numcells(:,1)) %if tetrode already listed in numcells
                    [row junk] = find(numcells(:,1) == tet);
                    numcells(row, 2) = numcells(row, 2) + length(cellinfo{d}{e}{tet});
                else
                    numcells = [numcells; tet length(cellinfo{d}{e}{tet}) ];
                end
            end
        end
    end
    [y i] = max(numcells(:,2));
    tetlist =  numcells(i,1);

elseif includeepochs == 2
    numcells = [0 0];
    d = day;
    for e = 1:length(task{d})
        if ~isempty(task{d}{e})
            if  isfield(task{d}{e},'type')
                if isequal(task{d}{e}.type,'run')
                    tetlist =  evaluatefilter(cellinfo{d}{e}, cellfilter);
                    % get rid of the cell indeces and extract only the tetrode numbers
                    if ~isempty(tetlist)
                        tetlist = unique(tetlist(:,1))';
                        for s = 1:length(tetlist)
                            tet = tetlist(s); %tetrode

                            if ismember(tet,numcells(:,1)) %if tetrode already listed in numcells
                                [row junk] = find(numcells(:,1) == tet);
                                numcells(row, 2) = numcells(row, 2) + length(cellinfo{d}{e}{tet});
                            else
                                numcells = [numcells; tet length(cellinfo{d}{e}{tet}) ];
                            end
                        end
                    end
                end
            end
        end
    end
    [y i] = max(numcells(:,2));
    tetlist =  numcells(i,1);

elseif includeepochs == 3
    numcells = [0 0];
    d = day;
    for e = 1:length(task{d})
        if ~isempty(task{d}{e})
            if  isfield(task{d}{e},'type')
                if isequal(task{d}{e}.type,'sleep')
                    tetlist =  evaluatefilter(cellinfo{d}{e}, cellfilter);
                    % get rid of the cell indeces and extract only the tetrode numbers
                    if ~isempty(tetlist)
                        tetlist = unique(tetlist(:,1))';
                        for s = 1:length(tetlist)
                            tet = tetlist(s); %tetrode

                            if ismember(tet,numcells(:,1)) %if tetrode already listed in numcells
                                [row junk] = find(numcells(:,1) == tet);
                                numcells(row, 2) = numcells(row, 2) + length(cellinfo{d}{e}{tet});
                            else
                                numcells = [numcells; tet length(cellinfo{d}{e}{tet}) ];
                            end
                        end
                    end
                end
            end
        end
    end
    [y i] = max(numcells(:,2));
    tetlist = numcells(i,1);
end

out = tetlist;
end