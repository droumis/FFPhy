function sj_matclust_isiviol(clustnum)

% 1) Look at ExcludeOvelappingPoints in matclust/ClustTools
% Can similarly exclude ISI Violations

%2) Option 2: Can filter out ISI VIolations 
% Ask for clustnum and isi in dialog. eg. ISIfilter in matclust/Filters
% Also look at other filters like allsleeps/allruns/burstfilter 
% out = true(length(clustdata.params(:,1)),1);
% isiviol_idx = find(diff(clustdata.params(clustattrib.clusters{clustnum}.index,1)) <=isi*10);
% isiviol_idx is the first of the 2 spikes. You can remove either isiviol_idx or isiviol_idx+1
% out(isiviol_idx+1) = false;


global clustdata;
global clustattrib;
global figattrib;

isi=2; % ms
isiviol_idx = find(diff(clustdata.params(clustattrib.clusters{clustnum}.index,1)) <=isi*10);

overlap = analyzeoverlap(0);
clustindex = find((overlap.list(:,1) == clustnum)|(overlap.list(:,2) == clustnum));

tmpadd = isiviol_idx;
tmpadd(:,2) = 0;
tmpadd = int32(tmpadd);
tmpadd(:,2) = fastbitset(tmpadd(:,2),2,true);
try
    clustattrib.eventeditindex{clustnum} = [clustattrib.eventeditindex{clustnum};tmpadd];
catch
    clustattrib.eventeditindex{clustnum} = tmpadd;
end

in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
in3 = false(length(clustdata.params(:,1)),1); %which points are excluded individually

try
    if ~isempty(clustattrib.eventeditindex{clustnum})
        in3(clustattrib.eventeditindex{clustnum}(:,1)) = true; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
    end
end

clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3));
matclust('addnewstate','exclude isi viol',figattrib.handles); 
matclust('plotgraph',figattrib.handles); 



%%% This is ExcludeOverlappingPoints
% overlap = analyzeoverlap(0);
% clustindex = find((overlap.list(:,1) == clustnum)|(overlap.list(:,2) == clustnum));
% for i = 1:length(clustindex)
%     currindex = overlap.points{i};
%     tmpadd = currindex(:);
%     tmpadd(:,2) = 0;
%     tmpadd = int32(tmpadd);
%     tmpadd(:,2) = fastbitset(tmpadd(:,2),2,true);
% 
%     try
%         clustattrib.eventeditindex{clustnum} = [clustattrib.eventeditindex{clustnum};tmpadd];
%     catch
%         clustattrib.eventeditindex{clustnum} = tmpadd;
%     end
%     
% end
% 
% in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
% in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
% in3 = false(length(clustdata.params(:,1)),1); %which points are excluded individually
% 
% try
%     if ~isempty(clustattrib.eventeditindex{clustnum})
%         in3(clustattrib.eventeditindex{clustnum}(:,1)) = true; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
%     end
% end
% 
% clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3));
% 
% matclust('addnewstate','exclude overlap',figattrib.handles); 
% matclust('plotgraph',figattrib.handles); 
