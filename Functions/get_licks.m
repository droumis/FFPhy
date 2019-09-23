

function F = get_licks(animal, index, DIO, task, varargin)

% create and save lick events;
saveout = 1;
if ~isempty(varargin)
    assign(varargin{:});
end

for x = 1:size(index,1)
    day = index(x,1);
    epoch = index(x,2);
    lickDIOIdx = task{day}{epoch}.inputdio;
    isinput = cellfun(@(x) isequal(x.input,1), DIO{day}{epoch}, 'un', 1);
    dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')),...
        DIO{day}{epoch}, 'un', 1);
    l = sortrows(cell2mat(cellfun(@(x) [x.times repmat(str2double(...
        regexp(x.original_id,'\d*','Match')), ...
        length(x.times),1)],DIO{day}{epoch}(ismember(dioID(isinput),lickDIOIdx)),'un',0)'),1); 
    lick{1,day}{1,epoch}.starttime = l(:,1);
    lick{1,day}{1,epoch}.id = l(:,2);
    lick{1,day}{1,epoch}.eventname = 'lick';
end
F = struct;
F(1).data = lick;
F(1).animal = animal;
if saveout
    save_data(F, 'filterframework', 'lick', 'animpos', 0, 'varname', 'lick');
end
end

    


