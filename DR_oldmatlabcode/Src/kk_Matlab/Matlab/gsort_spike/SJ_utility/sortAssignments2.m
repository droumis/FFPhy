function [reassignments, map] = sortAssignments2(assignments)

% SORTASSIGNMENTS  Renumbers assignments
%    reassignments = sortAssignments(assignments)
%
% Takes a list of assignment numbers and reassigns label numbers 
% such that the largest size group is assigned label '1', the next largest
% is assigned label '2', and so on.
% GET RID OF SIZE CONDITION IN THIS ONE:

clusters = unique(assignments);  % get a list of unique labels . . .
clusters(find(clusters == 0)) = []; % DO NOT WANT TO REASSIGN OUTLIERS BASED ON SIZE
numclusts = length(clusters);    %
clustsize = zeros(numclusts,1);  %

% SIZE OF CLUSTERS
for clust = 1:numclusts          % ... and count # elements assigned to each label
    clustsize(clust) = length(find(assignments == clusters(clust)));
end

% SORT BASED ON SIZE
% create a matrix with cols [old_label  num_elements] and (descending) sort on num_elemebts
 reassign_list = flipud(sortrows([clusters, clustsize], 2));

%%%%%%%%%%%% DEBUGGING - random assignments instead of size sorted.  Useful because
%%%%%%%%%%%%             it still gets rid of unused cluster numbers.
% reassign_list(1:numclusts,1) = reassign_list(randperm(numclusts),1);

% . . . and use that table to translate the original assignment list
reassignments = zeros(size(assignments));
for clust = 1:numclusts
    reassignments(find(assignments == clusters(clust,1))) = clust;
end
map = [unique(assignments) unique(reassignments)]; 

return;
