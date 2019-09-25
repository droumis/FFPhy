function f = gettetrodelist(f, cellfilter, includeepochs, varargin)
% tetlist = gettetrodelist(f,  cellfilter, includeepochs, 'options')
%
% takes filter f and changes f(an).data field from [tet cell] to [tet], selects
% tet witht he most cells using the following options
%
% cellfilter: example: '(isequal($area, ''CA3'')) '
% includeepochs, 1, 2, 3, default 1
%       1 -- all epochs
%       2 -- all run
%       3 -- all sleeps
% options:
%   ndays -- number of days cells were recorded on a particular
%           tetrode, default 1
%   maxcell -- 0, 1 ,2, default 0
%       0 include tetrodes regardless of number of cells
%       1  select a single tetrode for each day, select tetrode with most
%       cells recorded on that day
%       2  same as 1 but if all but 1 day use same tetrode, switch that day
%       to match all other days
%

ndays = 1;
maxcell = 0;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'ndays'
                ndays = varargin{option+1};
            case 'maxcell'
                maxcell = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

if maxcell == 0
    for an = 1:length(f)
        tetindex = [];
        for g = 1:length(f(an).data) %each epoch filter
            for e = 1:size(f(an).data{g},2) %each epoch

                tets = unique(f(an).data{g}{e}(:,1)); %find all tets with cells on them
                d = f(an).epochs{g}(e,1);
                tetindex = [tetindex; d*ones(length(tets),1) tets]; %day list of tetrodes
            end
        end
        tetindex = unique(tetindex, 'rows');
        table = tabulate(tetindex(:,2));
        count = table(:,1:2);
        tetlist{an} = count(count(:,2) >= ndays,1);
        for g = 1:length(f(an).data); %each epoch filter
            for e = 1:size(f(an).data{g},2) ;%each epoch
                f(an).data{g}{e} = tetlist{an};
            end
        end
    end


    %find tet with most cells
elseif maxcell ~= 0
    for an = 1:length(f)
        allmaxtet = [];
        cellinfo = loaddatastruct(f(an).animal{2},  f(an).animal{3}, 'cellinfo');
        task = loaddatastruct(f(an).animal{2},  f(an).animal{3},  'task');
        alldays = [];
        for g = 1:length(f(an).data) %each epoch filter
            alldays = [alldays; f(an).epochs{g}(:,1)];

        end
        alldays = unique(alldays);
        for i = 1:length(alldays)
            d = alldays(i);
            tet = gettetmaxcell(cellinfo, task, d, cellfilter, includeepochs);
            allmaxtet = [allmaxtet; d tet];
        end

        if maxcell == 2
            if length(unique(allmaxtet(:,2))) == 2;
                table = tabulate(allmaxtet(:,2));
                change = find(table(:,2)==1); %find tetrode that should be changed
                new = find(table(:,2)>1);
                allmaxtet(allmaxtet(:,2)==change,2)=new;
            end
        end
        tetlist{an} = allmaxtet;

        %change filter to reflect new tetrode assignments
        for g = 1:length(f(an).data) %each epoch filter
            for e = 1:size(f(an).data{g},2) %each epoch
                d= f(an).epochs{g}(e,1); %day
                tetrode = tetlist{an}(tetlist{an}(:,1)==d,2) ;%find the tetrode corresponding to day d
                f(an).data{g}{e} = tetrode;
            end
        end

    end
end
end



