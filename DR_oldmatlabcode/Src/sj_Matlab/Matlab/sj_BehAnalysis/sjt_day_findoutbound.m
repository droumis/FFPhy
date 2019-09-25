
function [outbound_logic, outbound_wellstend] = sjt_day_findoutbound (linpos,day,epoch)

% From linpos, finds start and end wells and whether correct and incorrect
% for all outbounbd trajectories on current day and epoch
% Shantanu Jadhav, 12/30/09

%
% Get Well Start and End for each time point
well_stend=linpos{day}{epoch}.statematrix.wellExitEnter;

% Get Trajectory start Indexes
aexit = find(diff(well_stend(:,1))~=0);
aenter = find(diff(well_stend(:,2))~=0);
trajst_idx = unique([aexit; aenter]);

% Get wells for unique trajectories
wells = well_stend(trajst_idx,:);

% Find outbound
outbound_stidx = find(wells(:,1)==1);
outbound_wellstend = wells(outbound_stidx,:);

% Find out which outbound is correct, and return the logic 
outbound_logic = zeros(length(outbound_wellstend),1);

mintr = min(outbound_stidx);
if mintr==1,
    outbound_stidx(find(outbound_stidx==1))=[];
end

corridx = find( (wells(outbound_stidx,2)~=wells(outbound_stidx-1,1)) & (wells(outbound_stidx,2)~=1) );
correct = outbound_stidx(corridx);


%% Either get rid of first trial, or keep it in as correct
% Alt 1 - Keep it as correct
% if mintr==1.
%     correct = [1;correct]; outbound_stidx = [1;outbound_stidx];
% end
% Alt 2 - Remove it - Do nothing and update length of outbound_logic
outbound_logic = zeros(length(outbound_wellstend)-1,1);

[c,ia,ib] = intersect(outbound_stidx, correct);

outbound_logic(ia) = 1;
sum(outbound_logic);



%%

% %% LONG WAY 
% 
% %% Get Well Start and End for each time point
% well_stend=linpos{day}{epoch}.statematrix.wellExitEnter;
% %% Find all times on Inbound trajectories 
% out_cum = find((well_stend(:,1)==1));
% %% Find Start Idxs of inbound trajectories
% d1=find(diff(out_cum)~=1);
% d2 = find( diff(well_stend(out_cum,2))~=0);
% d = unique([d1; d2]);
% outbound_stidx=([1; out_cum(d+1)]);
% 
% %% Find start and end wells of inbound trajectories
% outbound_wellstend = well_stend(outbound_stidx,:);
% returnwell_idx = find(outbound_wellstend(:,2)==1);
% %% Find out which inbound is correct, and return the logic in order of
% %% trial number
% outbound_logic = zeros(length(outbound_wellstend),1);
% corr = find(outbound_wellstend(:,2)==1);
% outbound_logic(corr) = 1;
