
function kk_updatetaskstruct_sleep(animdirect,prefix,days)

% Takes empty epochs in taskstruct and gives them type 'sleep'.
% Note that total epochs mirrors that of pos file.

for d=days
    
    taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix,d);
    load(taskfile);
    posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix,d);
    load(posfile);
    
    for e=1:length(pos{d})
        try
            if isempty(task{d}{e})
                task{d}{e}.type='sleep';
                disp(sprintf('%d %d converted to sleep',[d e]))
            end
        catch   % errors out past the last 'run' epoch -- assign the bookending sleep epoch
            task{d}{e}.type='sleep';
            disp(sprintf('%d %d converted to sleep',[d e]))
        end
    end
    
    
    save(taskfile,'task');
    
end


