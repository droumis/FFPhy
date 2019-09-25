function out = JY_calcinfieldtime(index, excludetimes, spikes, linpos, includestates, minV, varargin)
%
% returns the time indices the rat spends inside the placefield of each
% cell based on DFL_plotlinplacefielddata.m.
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
% OUTPUT:
%   []

close all;

day = index(1);
epoch_num = index(2);
tetrode_num = index(3);
neuron_num = index(4);

spatial_bins_per_seg = 50;

disp(sprintf('epoch %d, tetrode %d, neuron %d',epoch_num,tetrode_num,neuron_num));

% Get exclude indices
exclind = isExcluded(linpos{index(1,1)}{index(1,2)}.statematrix.time,excludetimes);

%Neural and position data to process
spike_data = spikes{index(1)}{index(2)}{index(3)}{index(4)};
linpos_data = linpos{index(1)}{index(2)};

timestep = linpos_data.statematrix.time(2) - linpos_data.statematrix.time(1);

% In linpos_data, set all segments within exclude times to be -1
linpos_data.statematrix.segmentIndex(find(~exclind)) = -1;


% bail out of function if spike_data.data is empty after filtering (this
% means there are no valid spikes to compute place fields for in this
% neuron)

if(isempty(spike_data.data))
    out = -1;
    return;
end

% Remove spikes that are in exclude times
spike_data.data(ismember(spike_data.data(:,7), find(~exclind)),:) = [];

%disp(sprintf('spikes excluded %d-%d=%d', length(spike_data.data),length(spike_data_data), length(spike_data.data) - length(spike_data_data) ));

% bail out of function if spike_data.data is empty after filtering (this
% means there are no valid spikes to compute place fields for in this
% neuron)

if(isempty(spike_data.data))
    out = -1;
    return;
end

% Occupancy table - one table for each epoch
occupancy_time = cell(1,length(linpos{1}));

total_neuron_count = 0;

normalized_firing_rate = cell(1,length(spikes{index(1)}{index(2)}));

% Setup occupancy table, array of segments (both directions) and bins per
% segment
occupancy_time = zeros(size(linpos_data.segmentInfo.segmentCoords,1),spatial_bins_per_seg);
temp_count = 0;
for seg_num = 1:size(occupancy_time,1)
    % Get indices of data on segment
    seg_data_ind = find(linpos_data.statematrix.segmentIndex == seg_num);
    temp_count = length(seg_data_ind) + temp_count;
    
    bin_length = linpos_data.segmentInfo.segmentLength(seg_num)/spatial_bins_per_seg;
    
    bin_pnt_count_total = 0;
    for bin_num = 1:size(occupancy_time,2)
        bin_coord = [(bin_num-1)*bin_length bin_num*bin_length];
        
        %count number of pos data points that falls within the bin_coord on this
        %segment
        lindist_seg_data = linpos_data.statematrix.linearDistanceToWells(seg_data_ind);
        bin_pnt_count = length(find(lindist_seg_data >= bin_coord(1) & lindist_seg_data < bin_coord(2)));
        
        bin_pnt_count_total = bin_pnt_count_total + bin_pnt_count;
        
        occupancy_time(seg_num, bin_num) = bin_pnt_count * timestep;
    end
end

seg_list = linpos_data.statematrix.segmentIndex(spike_data.data(:,7));


% Calculate normalized firing rate per segment bin
normalized_firing_rate = cell(1,length(occupancy_time));

spike_total = 0;
pos_total = 0;

% times points within placefields
timeindices=[];
% boundries of placefields
includeindices=[];

for seg_num = 1:size(occupancy_time,1)
    % Mask for the given segment applied to the pos/seg data
    seg_data_mask = linpos_data.statematrix.segmentIndex == seg_num;
    % Get indexes of data on segment
    bin_length = linpos_data.segmentInfo.segmentLength(seg_num)/spatial_bins_per_seg;
    
    % Loop through the same bins used in occupancy time
    norm_seg_spikes = zeros(1,size(occupancy_time,2));
    spike_seg_total = 0;
    
    % matrix to hold time indices

     
    for bin_num = 1:size(occupancy_time,2)
        % Calculate the coordinates of each bin
        bin_coord = [(bin_num-1)*bin_length bin_num*bin_length];
        % Apply segment mask to extract only segments we are current
        % processing in this loop
        lindist_seg_data = linpos_data.statematrix.linearDistanceToWells .* seg_data_mask;
        
        figure;
        list=find(lindist_seg_data~=0);
        plot(data{1,10}{1,2}.Pos.correcteddata(list,2),data{1,10}{1,2}.Pos.correcteddata(list,3),'.');
        
        
        % Narrow down points to bins
        pos_bin = find(lindist_seg_data > bin_coord(1) & lindist_seg_data <= bin_coord(2));
        
        % get the time indices for the positions
        currtimeindex=[];
        currtimeindex(:,1)=pos_bin;
        currtimeindex(:,2)=seg_num;
        currtimeindex(:,3)=bin_num;
        
        timeindices=[timeindices;currtimeindex];
        
        % find spikes that occurred on this point
        spike_count = length(find(ismember(spike_data.data(:,7),pos_bin)));
        
        % normalize spike count with occupancy time
        norm_seg_spikes(bin_num) = spike_count/occupancy_time(seg_num,bin_num);
        if(isnan(norm_seg_spikes(bin_num)))
            norm_seg_spikes(bin_num) = 0;
        end
        
        pos_total = pos_total + length(pos_bin);
        spike_seg_total = spike_seg_total + spike_count;
        
    end
    
    
    bin_coord = linspace(0,linpos_data.segmentInfo.segmentLength(seg_num)+linpos_data.segmentInfo.segmentLength(seg_num)/spatial_bins_per_seg,spatial_bins_per_seg);
    
    
    % gaussian smooth normalised spikes
    
    
    sigma=1;
    width = round(6./(linpos_data.segmentInfo.segmentLength(seg_num)/bin_num));
    support = (-width:width);
    gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
    gaussFilter = gaussFilter/ sum(gaussFilter);
    
    norm_seg_spikes_gauss=conv(norm_seg_spikes,gaussFilter,'same');
    
    normalized_firing_rate{seg_num} = norm_seg_spikes_gauss;
    
    % find bins with firing rate above 3Hz
    
    firingratethreshold=3;
    
    includebins=find(norm_seg_spikes_gauss > firingratethreshold);
    
    tempinclude=[ones(length(includebins),1)*seg_num' includebins'];
    includeindices=[includeindices;tempinclude];

end

% get all times when rat is in the placefield

inplacefieldtime=timeindices(ismember(timeindices(:,2:3), includeindices(:,1:2), 'rows'),1:2);



out.normalized_firing=normalized_firing_rate;
out.inplacefieldtimeindices=inplacefieldtime;






