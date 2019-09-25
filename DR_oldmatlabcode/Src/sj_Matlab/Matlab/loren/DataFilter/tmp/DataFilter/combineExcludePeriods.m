function out = combineExcludePeriods(excludeA, excludeB)
% out = combineExcludePeriods(excludeA, excludeB)
% Combines the exclude start and end times from excludeA and excludeB.
% Each input matrix is N by 2, where the colomns are the start and end
% times for each exclusion period.

tmpcombine = [excludeA; excludeB];
tmpborder = exclude2border(tmpcombine);
out = border2exclude(tmpborder);
%------------------------------------------------------
function out = border2exclude(borderMatrix)

statusVector = cumsum(borderMatrix(:,2));
statusVector(find(statusVector < 0)) = 0;

startperiods = borderMatrix((find(diff(statusVector) == -1)+1),1);
endperiods = borderMatrix((find(diff(statusVector) == 1)+1),1);

out = ([startperiods endperiods]); 
%------------------------------------------------------
function out = exclude2border(excludePeriods)


offVector = excludePeriods(:,1);
offVector(:,2) = -1;
onVector = excludePeriods(:,2);
onVector(:,2) = 1;

out = [[-inf 1]; sortrows([offVector;onVector],1)];