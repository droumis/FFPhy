
function [task] = sj_updatetaskenv(animdirect,prefix,day,epochs,env)

% Shantanu 
% To put env = 'env' in task struct, eg. 'lin', 'wtr1', 'wtr2'
% In addition to type such as 'sleep' in task struct which only has run types

% Task file
taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix, day);
load(taskfile);

for ep = epochs
    task{day}{ep}.environment = env;
end

% Save updated linpos file
save(taskfile,'task');
