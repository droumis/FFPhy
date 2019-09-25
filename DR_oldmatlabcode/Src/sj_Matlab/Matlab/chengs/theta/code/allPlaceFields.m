function pf= allPlaceFields(selectid)
%function pf= allPlaceFields(selectid)
%
% Return the locations of place fields.


setRoot
load([root '/data/analist-' selectid])
nC= length(analist.rat);
oldrat= ''; oldd= -1;
pf= nan*ones(nC, 2);
for iC= 1:nC
    % set up data
    rat= analist.rat{iC};
    num= analist.cellnum(iC,:); d=num(1); e=num(2); tet=num(3); c=num(4);
    if ~strcmp(oldrat, rat) 
        oldrat= rat; oldd= d;
        load(sprintf('/home/chengs/theta/%s/data2/select-%s.mat', rat, selectid));
    end
    nsel= getSelectId([d,e,tet,c], select);
    pf(iC,:)= select.a{nsel}{analist.iana(iC)}.linpos';
end
