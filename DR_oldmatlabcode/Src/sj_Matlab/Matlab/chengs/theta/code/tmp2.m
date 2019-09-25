function tmp2

global fmaux task spikedata
loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);

load(fmaux.select);
nsel= length(select.a);

for is=1:nsel
    num= select.cellnum(is,:);
    d=num(1); e=num(2); 
    arms= task{d}{e}.task;
    arms= arms(~ismember(arms, [1 3 7]));
    if isempty(arms)
        select.newarm(is,:)= nan;
    else
        select.newarm(is,:)= arms;
    end
    select.day(is,:)= task{d}{e}.exposure;
end

save([fmaux.select '2'], 'select');
