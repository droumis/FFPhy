function [newselect, totalana]= cleanSelect(select, confirm)
%function [newselect, totalana]= cleanSelect(select, confirm)
%  Remove empty entries from select-list
%
%  select       list to be cleaned
%  confirm= 0;  whether to manually confirm selection by visual inspection
%  newselect    cleaned list


if ~isfield(select, 'a') | ~isfield(select, 'cellnum')
    error('select must have fields ''cellnum'' and ''a'' ');
end

if nargin < 2; confirm= 0; end

newselect= [];
newn= 0;
totalana= 0;
for n=1:size(select.cellnum,1)

    if n > length(select.a) break; end
    if isempty(select.a{n}) continue; end

    saveCell= 1;

    if confirm
        d= select.cellnum(n,1);
        e= select.cellnum(n,2);
        t= select.cellnum(n,3);
        c= select.cellnum(n,4);

        showCellHist([d e t c]);
        in= input('save (y/n)?','s');
        if in== 'n'
            saveCell= 0;
        else 
            saveCell =1;
        end
    end
    if saveCell
        totalana= totalana+length(select.a{n});
        newn= newn+1;
        newselect.cellnum(newn,:)= select.cellnum(n,:);
        newselect.day(newn,:)= select.day(n,:);
        newselect.newarm(newn,:)= select.newarm(n,:);
        newselect.a{newn}= select.a{n};
        if isfield(select,'x') & n <= length(select.x)
            newselect.x{newn}= select.x{n};
        end
    end
end

