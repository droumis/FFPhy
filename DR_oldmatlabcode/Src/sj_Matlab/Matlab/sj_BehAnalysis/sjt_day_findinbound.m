
function [inbound_logic, inbound_wellstend] = sjt_day_findinbound (linpos,day,epoch)

% From linpos, finds start and end wells and whether correct and incorrect
% for all inbound trajectories on current day and epoch
% Shantanu Jadhav, 12/30/09


% Get Well Start and End for each time point
well_stend=linpos{day}{epoch}.statematrix.wellExitEnter;

% Get Trajectory start Indexes
aexit = find(diff(well_stend(:,1))~=0);
aenter = find(diff(well_stend(:,2))~=0);
trajst_idx = unique([aexit; aenter]);

% Get wells for unique trajectories
wells = well_stend(trajst_idx,:);

% Find inbound
inbound_stidx = find(wells(:,1)~=1);
inbound_wellstend = wells(inbound_stidx,:);
% Find out which inbound is correct, and return the logic 
inbound_logic = zeros(length(inbound_wellstend),1);
correct = find(inbound_wellstend(:,2)==1);
inbound_logic(correct) = 1;


%%

%%% LONG WAY %%%%
% %% Get Well Start and End for each time point
% well_stend=linpos{day}{epoch}.statematrix.wellExitEnter;
% %% Find all times on Inbound trajectories
% in_cum = find((well_stend(:,1)==3) | (well_stend(:,1)==2) );
% %% Find Start Idxs of inbound trajectories
% d1 = find(diff(in_cum)~=1);
% d2 = find( diff(well_stend(in_cum,1))~=0 | diff(well_stend(in_cum,2))~=0 );
% d = unique([d1; d2]);
% inbound_stidx=(in_cum(d+1));
% %% Find start and end wells of inbound trajectories
% inbound_wellstend = well_stend(inbound_stidx,:);
% %% Find out which inbound is correct, and return the logic in order of
% %% trial number
% inbound_logic = zeros(length(inbound_wellstend),1);
% corr = find(inbound_wellstend(:,2)==1);
% inbound_logic(corr) = 1;
