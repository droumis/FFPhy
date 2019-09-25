function out = DFL_plotlinplacefielddata(index, excludetimes, spikes, linpos, includestates, minV, f, varargin)
%
% out = DFL_plotlinplacefielddata(index, excludetimes, spikes, linpos, includestates, minV, f, varargin)
% Plot and saves the linearized place fields 
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

spatial_bins_per_seg = 20;

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
    % Get indexes of data on segment
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

% Plot histogram of number of spikes per segment
% subplot(7,2,1);
% set(gca,'FontSize',11);
% hist(seg_list,12)
% xlim([1 12]);
% title(sprintf('Epoch %d, Tetrode %d, Cell %d, spikes=%d',epoch_num, tetrode_num, neuron_num,size(spike_data.data,1)));
% xlabel('Segment #');
% ylabel('Spikes/Segment');

% Calculate normalized firing rate per segment bin
normalized_firing_rate_for = cell(1,length(occupancy_time));
normalized_firing_rate_rev = cell(1,length(occupancy_time));
spike_total_for = 0;
spike_total_rev = 0;
pos_total_for = 0;
pos_total_rev = 0;
for seg_num = 1:size(occupancy_time,1)
    % Mask for the given segment applied to the pos/seg data
    seg_data_mask = linpos_data.statematrix.segmentIndex == seg_num;
    % Get indexes of data on segment
    bin_length = linpos_data.segmentInfo.segmentLength(seg_num)/spatial_bins_per_seg;
      
    
    % get the direction of all the time points, direction is stored in
    % statematrix.wellExitEnter
    
    seg_datamaskindex=find(seg_data_mask==1);
    seg_directions=find(linpos_data.statematrix.wellExitEnter(seg_data_mask==1,:));
    uniquedirections=unique(seg_directions,'rows');
    directionindex1=seg_data_mask(rowfind(uniquedirections(1,:),seg_directions)==1);
    directionindex2=seg_data_mask(rowfind(uniquedirections(2,:),seg_directions)==1);
    
  
    % Loop through the same bins used in occupancy time
    norm_seg_spikes_for = zeros(1,size(occupancy_time,2));
    norm_seg_spikes_rev = zeros(1,size(occupancy_time,2));
    spike_seg_total = 0;
    for bin_num = 1:size(occupancy_time,2)
        % Calculate the coordinates of each bin
        bin_coord = [(bin_num-1)*bin_length bin_num*bin_length];
        % Apply segment mask to extract only segments we are current
        % processing in this loop
        lindist_seg_data_for = linpos_data.statematrix.linearDistanceToWells .* directionindex1;
        lindist_seg_data_rev = linpos_data.statematrix.linearDistanceToWells .* directionindex2;
        
        
        % Narrow down points to bins
        pos_bin_for = find(lindist_seg_data_for > bin_coord(1) & lindist_seg_data_for <= bin_coord(2));
        pos_bin_rev = find(lindist_seg_data_rev > bin_coord(1) & lindist_seg_data_rev <= bin_coord(2));
        
        
        % find spikes that occurred on this point
        spike_count_for = length(find(ismember(spike_data.data(:,7),pos_bin_for)));
        spike_count_rev = length(find(ismember(spike_data.data(:,7),pos_bin_rev)));
        
        % normalize spike count with occupancy time
        norm_seg_spikes_for(bin_num) = spike_count_for/occupancy_time(seg_num,bin_num);
        norm_seg_spikes_rev(bin_num) = spike_count_rev/occupancy_time(seg_num,bin_num);
        
        if(isnan(norm_seg_spikes_for(bin_num)))
            norm_seg_spikes_for(bin_num) = 0;
        end
        
        if(isnan(norm_seg_spikes_rev(bin_num)))
            norm_seg_spikes_rev(bin_num) = 0;
        end
        
        
        pos_total_for = pos_total_for + length(pos_bin_for);
        pos_total_rev = pos_total_rev + length(pos_bin_rev);
        
        spike_seg_total_for = spike_seg_total_for + spike_count_for;
         spike_seg_total_rev = spike_seg_total_rev + spike_count_rev;
    end
    bin_coord = linspace(0,linpos_data.segmentInfo.segmentLength(seg_num)+linpos_data.segmentInfo.segmentLength(seg_num)/spatial_bins_per_seg,spatial_bins_per_seg);
    normalized_firing_rate_for{seg_num} = norm_seg_spikes_for;
    normalized_firing_rate_rev{seg_num} = norm_seg_spikes_rev;
    
    % Plot firing rates per segment
    figure(1);
    set(gcf,'position',[0 0 800 900]);
    set(gcf,'PaperPositionMode','auto');
    subplot(7,2,seg_num+1)
    set(gca,'FontSize',11);
    plot(bin_coord,normalized_firing_rate_for{seg_num},'b');
    hold on;
    plot(bin_coord,normalized_firing_rate_rev{seg_num},'r');
    
    title(sprintf('Firing Rate on Segment (D%d E%d T%d N%d S%d n=%d)',day, epoch_num, tetrode_num, neuron_num, seg_num, spike_seg_total));
    xlabel('Segment Distance');
    ylabel('Firing Rate');
    
    
    % Visualize the firing rate using the special plot_segment_linpf
    % function
%     figure(2);
%     set(gca,'FontSize',14);
%     endpt1 = linpos_data.segmentInfo.segmentCoords(seg_num,1:2);
%     endpt2 = linpos_data.segmentInfo.segmentCoords(seg_num,3:4);
%     plot_segment_linpf(endpt1,endpt2,bin_coord,normalized_firing_rate{seg_num});
%     title(sprintf('Linearized Place Field (D%d E%d T%d N%d)', day, epoch_num, tetrode_num, neuron_num));
%     xlabel('x pos');
%     ylabel('y pos');
end

% loop through all plots and set ylim to same
max_firing_rate = max(cellfun(@(x) (max(x)), normalized_firing_rate));
if max_firing_rate == 0
    max_firing_rate = 0.01;
end

figure(1)
for seg_num = 1:size(occupancy_time,1)
    subplot(7,2,seg_num+1)
    ylim([0,max_firing_rate]);
end

figure(1)
[s,mess,messid] = mkdir(sprintf('%sPlot/linplacefields/',f.animal{2}));
print(sprintf('%sPlot/linplacefields/%s_linpf_d%de%dt%dc%d',f.animal{2},f.animal{1},index(1),index(2),index(3),index(4)),'-depsc');

% figure(2)
% print(sprintf('%sPlot/linplacefields/%s_linpf_d%de%dt%dc%d_maze',f.animal{2},f.animal{1},index(1),index(2),index(3),index(4)),'-depsc');


out = 1;





% Special plot function to use to plot linear place fields and visualize on
% maze
function plot_segment_linpf(a,b,bin_coord,normalized_firing_rate)
color_map = hsv(20);

if b(1)-a(1) < 0
    a_temp = a;
    a = b;
    b = a_temp;
    normalized_firing_rate = normalized_firing_rate(end:-1:1);
end

hold off;
line([a(1) b(1)],[a(2) b(2)]);
hold on;
pt_prev = a;
for ii = 2:length(bin_coord)
    % calculate endpoint coordinates of next segment
    pt_next = [bin_coord(ii-1)*cos(atan((b(2)-a(2))/(b(1)-a(1)))) + a(1), bin_coord(ii-1)*sin(atan((b(2)-a(2))/(b(1)-a(1)))) + a(2)];
    % Draw colored line segment to represent firing on that bin of the
    % segment
    line([pt_prev(1) pt_next(1)],[pt_prev(2) pt_next(2)],'Color',ind2rgb(floor((normalized_firing_rate(ii-1)))+1,jet(20)),'LineWidth',10);
    pt_prev = pt_next;
end
