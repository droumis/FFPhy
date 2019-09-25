
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


%dothis=1;
% Across all run epochs
for i = 1:length(task{day})
    if ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'run')) )
        
        disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
        
        % Check if linpos has already been updated for this day and epoch
        %if dothis==1
        if (~isfield(linpos{day}{i}.statematrix,'traj')) || (~isfield(linpos{day}{i}.statematrix,'lindist'))
            % If either field is missing, check if behavestate file exists,
            if exist(statefile,'file');
                % If yes, update linpos
                linpos{day}{i}.statematrix.traj = behavestate{day}{i}.state;
                linpos{day}{i}.statematrix.lindist = behavestate{day}{i}.lindist;
            else
                % if not, calculate behavestate and then update
                [state, lindist] = sj_getbehavestate(linpos, day, i, includeStates, 'minlinvelocity', minvelstate);
                linpos{day}{i}.statematrix.traj = state;
                linpos{day}{i}.statematrix.lindist = lindist;
                behavestate{day}{i}.state=state;
                behavestate{day}{i}.lindist=lindist;
            end
        end
    end
end

% Save updated linpos file
save(linposfile,'linpos');

% Save behavestate if it did not exist originally
if ~exist(statefile,'file');
    save(statefile,'behavestate');   
end

