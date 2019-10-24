

function F = get_licks(animal, index, DIO, task, varargin)

% create and save lick events;
% example of script that calls this: 
% lickXcorrSU_20190916.m
saveout = 1;
if ~isempty(varargin)
    assign(varargin{:});
end

for x = 1:size(index,1)
    day = index(x,1);
    epoch = index(x,2);
    % from lick DIO ID in taskInfo to index into DIO
    lickDIOID = task{day}{epoch}.inputdio;
    isinput = cellfun(@(x) isequal(x.input,1), DIO{day}{epoch}, 'un', 1);
    dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')),...
        DIO{day}{epoch}, 'un', 1);
    lickDIOIdx = find(ismember(dioID(isinput),lickDIOID));
    lickPortDin = DIO{day}{epoch}(lickDIOIdx);
    % for each lick port Din, gather it's times along with index
    l = cellfun(@(x) x.times, lickPortDin,'un',0)';
    [licktime,X] = sort(cell2mat(l));
    
    for d = 1:length(lickDIOID)
        g{d,1} = repmat(lickDIOID(d), length(l{d}),1);
    end
    i = cell2mat(g);
    id = i(X);
    lick{1,day}{1,epoch}.eventtime = licktime;
    lick{1,day}{1,epoch}.id = id;
    lick{1,day}{1,epoch}.eventname = 'lick';
end
F = struct;
F(1).data = lick;
F(1).animal = animal;
if saveout
    save_data(F, 'filterframework', 'lick', 'animpos', 0, 'varname', 'lick');
end
end

    


