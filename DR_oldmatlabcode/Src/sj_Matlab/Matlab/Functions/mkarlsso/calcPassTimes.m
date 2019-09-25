function out = calcPassTimes(index, excludetimes, linpos, varargin)

% calulates the amount of time it takes to get from one well to another,
% for each pass between two different wells.
statematrix = linpos{index(1)}{index(2)}.statematrix;

passThresh = 5;  %Linear distance from well that starts and ends a pass 
inPass = ( (statematrix.linearDistanceToWells(:,1) > passThresh) & ...
           (statematrix.linearDistanceToWells(:,2) > passThresh) & ...
           (statematrix.linearDistanceToWells(:,3) > passThresh) );
       
%find all the times when the trajectory changes
wellDiff = diff(statematrix.wellExitEnter);
trajChange = find((wellDiff(:,1) ~= 0) | (wellDiff(:,2) ~= 0));

passtimes = [];
%calculate the time spent for each pass
for i = 1:length(trajChange)-1
    startInd = trajChange(i)+1;
    endInd = trajChange(i+1);
    trajInPass = find(inPass(startInd:endInd));
    differentwells = (statematrix.wellExitEnter(startInd,1) ~= statematrix.wellExitEnter(startInd,2));
    if ( (~isempty(trajInPass)) & differentwells)
        newStartInd = trajInPass(1)+startInd-1;
        newEndInd = trajInPass(end)+startInd-1;
        passtimes = [passtimes; (statematrix.time(newEndInd) - statematrix.time(newStartInd))];
    end
end

out.passtimes = passtimes;

    