function [out] = dfa_getSpikePos(index, excludeperiods, spikes,linpos, pos, task, varargin)

if (~isempty(varargin))
    assign(varargin{:});
end

warning('OFF','MATLAB:divideByZero');
out.index = index;

try spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
catch
    out.spikes = [];
    return
end

spikesData = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
posdata = pos{index(1)}{index(2)}.data;
posfields = pos{index(1)}{index(2)}.fields;
taskEnv = task{index(1)}{index(2)}.environment;

timestring = 'time';
timecol = find(cell2mat(cellfun(@(x) strcmp(x,timestring), strsplit(posfields, ' '), 'UniformOutput', false)));
postime = posdata(:,timecol);

if strcmp(taskEnv, 'wtrack')
    statematrix = linpos{index(1)}{index(2)}.statematrix;
elseif strcmp(taskEnv, 'openfield')
    statematrix.time = postime;
    statematrix.traj = ones(length(postime),1);
end

statevector = statematrix.traj;

% front padding the posdata vect to match statevec
if length(statevector(:,1)) ~= length(posdata(:,1));
    posdata = [zeros((length(statevector(:,1))-length(posdata(:,1))), length(posdata(1,:))); posdata];
end

% Use Exclude Periods for TimeFilter version in addition to statevector=-1
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1; % Based on exclude time,

posindexfield = knnsearch(postime, spikesData(:,1));

xstring = 'x-sm';
xcol = find(cell2mat(cellfun(@(x) strcmp(x,xstring), strsplit(posfields, ' '), 'UniformOutput', false)));
posxdata = posdata(posindexfield,xcol);
ystring = 'y-sm';
ycol = find(cell2mat(cellfun(@(x) strcmp(x,ystring), strsplit(posfields, ' '), 'UniformOutput', false)));
posydata = posdata(posindexfield,ycol);

if ~isempty(spikesData)
    spikesData = [spikesData(:,1) posxdata posydata];
else
    spikesData = [];
end
out.spikes = spikesData;
out.spikesfields = {'spiketime', 'posx', 'posy'};
end