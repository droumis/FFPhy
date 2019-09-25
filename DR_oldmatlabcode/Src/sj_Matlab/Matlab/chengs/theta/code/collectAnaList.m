function analist= collectAnaList(selectid)
%function analist= collectAnaList(selectid)

prefix={'kyl', 'ter', 'sta', 'fel'};

pro_retro= 0;

global fmaux % file manager auxillary variable
olddir= pwd;
nRats= length(prefix);
analist= {};
if pro_retro
    analist.proind= [];
    analist.retroind= [];
end
nA= 0;
%analist.adaptid= adaptid;
analist.selectid= selectid;
analist.ncells= 0;
analist.rat={};
analist.cellnum=[];
analist.traj=[];
analist.iana=[];
for r=1:nRats
    % change directory, set local options
    cd(['/home/chengs/theta/' prefix{r} '/data2']);
    setLocalOptions;

    % select file
    fmaux.selectid= selectid;
    fmaux.select=[fmaux.data2dir '/select-' fmaux.selectid];
    load(fmaux.select);

    analist.(['ind' prefix{r} ])= [];
    if pro_retro
        analist.([prefix{r} 'proind'])= [];
        analist.([prefix{r} 'retroind'])= [];
        analist.ncells= analist.ncells + size(select.cellnum,1);
        for iC=1:size(select.cellnum,1)
            for j= 1:2:length(select.a{iC});
                nA= nA+1;
                analist.rat{nA,1}= prefix{r};
                analist.cellnum(nA,:)= select.cellnum(iC,:);
                analist.day(nA,:)= select.day(iC,:);
                analist.newarm(nA,:)= select.newarm(iC);
                analist.nsel(nA,:)= iC;
                analist.iana(nA,:)= j;
                analist.(['ind' prefix{r} ])(end+1)= nA;
                analist.traj(nA,:)= [select.a{iC}{j}.traj, select.a{iC}{j+1}.traj];
                if find(analist.traj(nA)==0)
                    analist.proind(end+1,1)= nA;
                    analist.([prefix{r} 'proind'])(end+1,1)= nA;
                else 
                    analist.retroind(end+1,1)= nA;
                    analist.([prefix{r} 'retroind'])(end+1,1)= nA;
                end
            end
        end
    else
        analist.ncells= analist.ncells + size(select.cellnum,1);
        for iC=1:size(select.cellnum,1)
            for j= 1:length(select.a{iC});
                nA= nA+1;
                analist.rat{nA,1}= prefix{r};
                analist.cellnum(nA,:)= select.cellnum(iC,:);
                analist.day(nA,:)= select.day(iC);
                analist.newarm(nA,:)= select.newarm(iC);
                analist.nsel(nA,:)= iC;
                analist.iana(nA,:)= j;
                analist.([ 'ind' prefix{r}])(end+1,1)= nA;
                analist.traj(nA,:)= select.a{iC}{j}.traj;
            end
        end
    end
%    keyboard
end
cd(olddir)

save(['/home/chengs/theta/data/analist-' selectid], 'analist');
%save(['/home/chengs/theta/data/analist-' selectid '-' adaptid], 'analist');


