[cl,nc]= findDiffCells('pf7-novelArm', 'pf6-novelArm');
load analist-pf7-novelArm
%[cl,nc]= findDiffCells('placefields4-novel2', 'pf7-novelArm');
%load analist-placefields4-novel2
%[cl,nc]= findDiffCells('placefields4-train2', 'pf7-fam');
%load analist-placefields4-train2

al=[];
na= length(analist.rat);
ja= 0;
for ia=1:na
    j= find(strcmp(cl.rat, analist.rat{ia})' &...
        cl.cellnum(:,1)==analist.cellnum(ia,1)&...
        cl.cellnum(:,2)==analist.cellnum(ia,2)&...
        cl.cellnum(:,3)==analist.cellnum(ia,3)&...
        cl.cellnum(:,4)==analist.cellnum(ia,4) );
    if ~isempty(j)
        ja=ja+1;
        al.rat{ja}= analist.rat{ia};
        al.cellnum(ja,:)= analist.cellnum(ia,:);
        al.traj(ja,:)= analist.traj(ia,:);
        al.newarm(ja,:)= analist.newarm(ia,:);
        al.day(ja,:)= analist.day(ia,:);
        al.iana(ja,:)= analist.iana(ia,:);
    end
end
al.ncells= nc;
analist= al;
save analist-novelArmDiff_7_6 analist
%save analist-placefields4-trainDiff analist
