%function anaMassive
%function anaMassive(selectid)

nspikes= 20;
dt= 0.25;
%selectid= 'CA1PE';
%selectid= 'pf7';
%selectid= 'pf7-fam';
%selectid= 'pf7-novelArm';
selectid= 'pf7-famArm';
%selectid= 'placefields4-fam2';
%selectid= 'placefields4-novel2';
%selectid= 'placefields4-train2';
[cl,nc]= collectCellList(selectid);

nevents= zeros(nc,1);
Levents= zeros(nc,1);
setRoot;
for ic=1:nc
    monitorProgress(ic,nc);
    rat= cl.rat{ic};
    d= cl.cellnum(ic,1); e= cl.cellnum(ic,2); 
    t= cl.cellnum(ic,3); c= cl.cellnum(ic,4);
    if ic==1 | ~strcmp(rat, cl.rat{ic-1}) | d~= cl.cellnum(ic-1,1)
        load(sprintf('%s/%s/data2/spikedata%.2d', root, rat, d));
    end
    st= spikedata{d}{e}{t}{c}.time;
    nsp= length(st);
    sep= st(nspikes:nsp)-st(1:nsp-nspikes+1);
    ind= find(sep<=0.5);
    [lo,hi]= findcontiguous(ind,nspikes);
    if isempty(lo); continue; end
    nevents(ic)= length(lo);
    Levents(ic)= mean(hi-lo);
end


for d=1:4
    if d<4
        ind= find(cl.day==d);
    else
        ind= find(cl.day>10);
    end
    N(d)= sum(nevents(ind));
    ncells(d)= length(ind);
    meanN(d)= N(d)/ ncells(d);
    stdN(d)= std(nevents(ind));
    seN(d)= stdN(d)/sqrt(ncells(d));
end

%N,ncells,meanN,seN
disp(selectid)
disp(meanN)
disp(seN)
