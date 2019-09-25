function myselect(inid, outid)
% function myselect(infile, outfile)

global fmaux
load([fmaux.data2dir '/select-' inid]);

newselect= {};
newn= 0;
for i=1:length(select.cellnum)
    good= 0;
    d= select.cellnum(i,1);
    e= select.cellnum(i,2);
    t= select.cellnum(i,3);
    c= select.cellnum(i,4);

    % begin selection criteria
    if(e== 4) good= 1; end

    % end selection criteria

    if(good)
        newn= newn+1;
        newselect.cellnum(newn,:)= select.cellnum(i,:);
        if(isfield(select,'a'))
            newselect.a{newn}= select.a{i};
        end
        if(isfield(select,'x'))
            newselect.x{newn}= select.x{i};
        end
    end
end
select= newselect;
save([fmaux.data2dir '/select-' inid], 'select');
