
function kk_updatetaskstruct_environment(animdirect,prefix,index,environment)

% assigns environment field to selected days and epochs
% the index matrix are day - epoch pairs (as in createtaskstruct)
% so day 3 epoch 2 and day 4 epoch 1 is [3 2 ; 4 1]

for p=1:size(index,1)
    day=index(p,1);
    epoch=index(p,2);
    taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix, day);
    load(taskfile);
    task{day}{epoch}.environment = environment;
    % Save updated linpos file
    save(taskfile,'task');
end

