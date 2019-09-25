%load select-all-placefields
load select-all-placefields-ProRetro
select
nana=length(select.a);
for n=1:nana
    fprintf(1,'%d: ---\n',n);
    nt= length(select.a{n});
    for i=1:nt
        fprintf(1, '   traj= %d',select.a{n}{i}.traj);
        disp(select.a{n}{i}.linpos');
    end
end
