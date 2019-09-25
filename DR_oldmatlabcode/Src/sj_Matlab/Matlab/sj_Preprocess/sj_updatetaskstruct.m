
function [task] = sj_updatetaskstruct(animdirect,prefix,day,epochs,type)

% Shantanu Jun2011
% To put type = 'sleep' in task struct which only has run types

% Task file
taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix, day);
load(taskfile);

for ep = epochs
    task{day}{ep}.type = type;
end

% Save updated linpos file
save(taskfile,'task');
