function [tetrodes, epochs]= collectTet(selectid)
%function [tetrodes, epochs]= collectTet(selectid)
%
% Assemble a list of all tetrodes and epochs.
%

global task tetloc
tetrodes= []; ntet= 0;
epochs= []; nepoch= 0;

if nargin<1; selectid= 'all'; end
%    load(['/home/chengs/theta/data/analist-' selectid]);
%    al= analist;
%    nc= length(al.rat);
[al,nc]= collectCellList(selectid);
for ic=1:nc
    rat= al.rat{ic};
    d= al.cellnum(ic,1); e= al.cellnum(ic,2); t= al.cellnum(ic,3);
    if(ic==1 | ~strcmp(rat, al.rat{ic-1}) | d~=al.cellnum(ic-1,1)...
        | e~=al.cellnum(ic-1,2) | t~=al.cellnum(ic-1,3))
        ntet= ntet+1;
        tetrodes.num(ntet, :)= [d e t];
        tetrodes.day(ntet,1)= al.day(ic);
        tetrodes.newarm(ntet,1)= al.newarm(ic);
        tetrodes.rat{ntet}= al.rat{ic};
        tetrodes.ncells(ntet,:)= sum(strcmp(rat, al.rat)' & d==al.cellnum(:,1) & ...
            e==al.cellnum(:,2) & t==al.cellnum(:,3));
    end

    if(ic==1 | ~strcmp(rat, al.rat{ic-1}) | d~=al.cellnum(ic-1,1)...
        | e~=al.cellnum(ic-1,2))
        nepoch= nepoch+1;
        epochs.num(nepoch, :)= [d e];
        epochs.day(nepoch, :)= al.day(ic);
        epochs.newarm(nepoch, :)= al.newarm(ic);
        epochs.rat{nepoch}= al.rat{ic};
        epochs.ncells(nepoch,:)= sum(strcmp(rat, al.rat)' & d==al.cellnum(:,1) & ...
            e==al.cellnum(:,2));
    end
end
save(['/home/chengs/theta/data/tetrodes-' selectid], 'tetrodes');
save(['/home/chengs/theta/data/epochs-' selectid], 'epochs');


