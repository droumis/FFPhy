function n= getAnaListId(selectid, rat, num)
%function n= getAnaListId(selectid, rat, num)
% n can be a vector!

setRoot;
load([root '/data/analist-' selectid])
nC= length(analist.rat);
n= [];
nfound= 0;
for iC= 1:nC
    % set up data
    if strcmp(rat, analist.rat{iC});
        if all(num== analist.cellnum(iC,:))
            nfound= nfound+1;
            n(nfound,1)= iC;
        end
    end
end

