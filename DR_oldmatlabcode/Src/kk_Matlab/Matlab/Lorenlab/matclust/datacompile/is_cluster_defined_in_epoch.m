function out = is_cluster_defined_in_epoch(cluster,epoch)

%IS_CLUSTER_DEFINED_IN_EPOCH(CLUSTER,EPOCH)
%Assumes that a MATCLUST file is currently open.
%This function returns a 1 if the cluster has been defined for any
%points within the time epoch.  If not, it returns a zero.  If the
%specified time epoch or cluster does not exist in the matclust file, it returns an
%empty array.  Note: epoch 1 is always the '"Allpoints" time filter, and epoch 2
%is the first user-defined time filter, and so on.  This function works
%correctly even when the user has created more complicated filters.  For
%example, if the user has created an "other" filter that includes only points
%from epoch 2 and epoch 5, and decided to cluster something only in this
%filter, then this program will recognize that that the cluster was defined
%for epochs 2 and 5, but not 3 or 4.  (The cluster will be of course be defined
%for epoch 1 because this includes all points).

global clustattrib
global clustdata

defallpoints = 0;
defepoch = 0;

defineaxes = clustattrib.clusters{cluster}.defineaxes;

try
    definedfilters = defineaxes(:,3);
    epochtimes = clustdata.timefilterranges(epoch,1:2);
catch
    out = [];  %if epoch if not defined, return an empty array
    return
end
for i = 1:length(definedfilters)
    tmptimeepoch = clustattrib.filterindex{definedfilters(i)}(1);
    tmpotherfilts = clustattrib.filterindex{definedfilters(i)}(2:end);
    tmpfilter = fastandbit(clustdata.timefilters,tmptimeepoch);
    tmpfilter2 = fastandbit(clustdata.otherfilters,tmpotherfilts);

    tmppoints = find(tmpfilter & tmpfilter2);
    times = clustdata.params(tmppoints,1);
    pointsinepoch = find((clustdata.params(:,1) >= epochtimes(1))&(clustdata.params(:,1) < epochtimes(2)));
    if ~isempty(intersect(pointsinepoch,tmppoints))
        defepoch = 1;
        break;
    end
    

%      if (tmptimeepoch == 1)
%          defallpoints = 1;
%      end
%      if (tmptimeepoch == epoch)
%          defepoch = 1;
%      end
end

if (defallpoints|defepoch)
    out = 1;
else
    out = 0;
end

