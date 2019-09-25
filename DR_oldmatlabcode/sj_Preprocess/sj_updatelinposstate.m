
function [linpos] = sj_updatelinposstate(animdirect,prefix,day)

% Shantanu Jun2011
% If linpos had not been updated with getbehavestate to get statematrix.traj and statematrix.lindist, do it now. 
% For old versions of linpos files

% For getbehavestate
includeStates = 6;
minvelstate = 3; %cm/sec

% Linpos file
linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
load(linposfile);
% Task file
taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix, day);
load(taskfile);
% Load behavestate file if exists
statefile = sprintf('%s/%sbehavestate%02d.mat', animdirect, prefix, day);
if exist(statefile,'file');
    load(statefile);   
end


% Across all run epochs
for i = 1:length(task{day})   
      if ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'run')) )
          
          disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
          
          if exist(statefile,'file');
              linpos{day}{i}.statematrix.traj = behavestate{day}{i}.state;
              linpos{day}{i}.statematrix.lindist = behavestate{day}{i}.lindist;
          else
              [state, lindist] = sj_getbehavestate(linpos, day, i, includeStates, 'minlinvelocity', minvelstate);
              linpos{day}{i}.statematrix.traj = state;
              linpos{day}{i}.statematrix.lindist = lindist;
              behavestate{day}{i}.state=state;
              behavestate{day}{i}.lindist=lindist;  
          end
      end     
end

% Save updated linpos file
save(linposfile,'linpos');

% Save behavestate if it did not exist originally
if ~exist(statefile,'file');
    save(statefile,'behavestate');   
end

