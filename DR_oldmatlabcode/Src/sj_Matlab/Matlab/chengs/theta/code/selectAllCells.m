function selectAllCells
%function selectAllCells
%  
%  Build database of all cells with at least one recorded spike.

global fmaux task spikedata
loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);

ndays = length(task);
nsel=0;
select=zeros(1,4); 

% select all cells with at least one recorded spike
for d = 1:ndays
    if isempty(task{d}); continue; end
    loadVar(fmaux.data2dir, 'spikedata', d);
    for e = 1:length(task{d})
        if isempty(task{d}{e}) | ~strncmp(task{d}{e}.type, 'run', 3);
            continue;
        end
        for t = 1:length(spikedata{d}{e})
            for k = 1:length(spikedata{d}{e}{t}) 
                if (~isempty(spikedata{d}{e}{t}{k})) &...
                    ~isempty(spikedata{d}{e}{t}{k}.time)
                    nsel=nsel+1;
                    select.cellnum(nsel,:)=[d e t k];
                    select.day(nsel,:)= task{d}{e}.exposure;
                    arms= task{d}{e}.task;
                    arms= arms(~ismember(arms, [1 3 7]));
                    if isempty(arms)
                        select.newarm(nsel,:)= nan;
                    else
                        select.newarm(nsel,:)= arms;
                    end
                    for j=1:4
                        select.a{nsel}{j}.traj= j-1;
                        select.x{nsel}{j}= {};
                    end
                end
            end
        end
    end
end

if nsel==0; select=[]; end

save([fmaux.data2dir '/select-all'], 'select');
