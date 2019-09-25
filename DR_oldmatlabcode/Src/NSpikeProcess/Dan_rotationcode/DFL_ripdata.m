function [out] = DFL_ripdata(index, excludetimes, spikes, data, includestates, minV, f, varargin)
% out = DFL_ripdata(index, excludetimes, spikes, includestates, minV, f, varargin)
% Filter function that counts the number of unique cells firing in a given
% window and saves it to a file.  Also includes a function that will try to
% display the cells firing within each ripple.
% Recommended Filter: multicellanal
% Filter Inputs:
%   spikes - spike structure from filter framework
% Options:
%   None
% OUTPUT:
%   []

%close all;

% Assume day and epoch are all the same
day = index(1,1);
epoch = index(1,2);
dsz = '';
if (day < 10)
    dsz = '0';
end
filename=strcat(f.animal{2},f.animal{1},'ripple_cell_count',dsz,num2str(day),'.mat');
if(exist(filename,'file'));
    load(filename);
end

ripple_spikes = cell(1,day);
ripple_spikes{1,day} = cell(1,epoch);
ripple_spikes{1,day}{1,epoch} = cell(2,size(excludetimes,1));

%iterates through all cells
for ii = 1:size(index,1)
    
    tetrode_num = index(ii,3);
    neuron_num = index(ii,4);
    
    
    disp(sprintf('epoch %d, tetrode %d, neuron %d',epoch,tetrode_num,neuron_num));
    
    
    %Neural and position data to process
    spike_data = spikes{day}{epoch}{tetrode_num}{neuron_num};
    
    %timestep = linpos_data.statematrix.time(2) - linpos_data.statematrix.time(1);
    
    % get all spikes within actual epoch according to position
    % reconstruction
    epochst=data{day}{epoch}.Pos.correcteddata(1,1);
    epochen=data{day}{epoch}.Pos.correcteddata(end,1);
    
    spike_data.data=spike_data.data(find(isExcluded(spike_data.data(:,1),[epochst epochen])),:);
    
    
    
    
    % bail out of iteration if spike_data.data is empty after filtering (this
    % means there are no valid spikes to compute place fields for in this
    % neuron)
    
    if(isempty(spike_data.data))
        out = -1;
        continue;
    end
    
    % bail out of iteration if spike_data.data is empty after filtering (this
    % means there are no valid spikes to compute place fields for in this
    % neuron)
    if(isempty(spike_data.data))
        out = -1;
        continue;
    end
    
    
    
    for ind = 1:size(excludetimes,1)-1
        rippleind = find(spike_data.data(:,1) > excludetimes(ind,2) & spike_data.data(:,1) < excludetimes(ind+1,1));
        % start time of ripple
        ripple_spikes{day}{epoch}{1,ind} = excludetimes(ind,2);
        if(~isempty(rippleind))
            ripple_spikes{day}{epoch}{2,ind} = [ripple_spikes{day}{epoch}{2,ind};spike_data.data(rippleind,:), repmat([tetrode_num, neuron_num],length(rippleind),1)];
        end
    end
    
    
end

% is the barrier up for the ripple
if ~isempty(data{day}{epoch}.Events.Barrier)
barriertime=[data{day}{epoch}.Events.Barrier(1,1) data{day}{epoch}.Events.Barrier(1,1)+data{day}{epoch}.Events.Barrier(1,4)]./10000;
ripple_barrier=isExcluded(excludetimes(:,1),barriertime);
else
    ripple_barrier=zeros(size(excludetimes,1),1);
end

% sort spikes in ripples by their start times
for ind = 1:size(excludetimes,1)-1
    ripple_spikes{day}{epoch}{2,ind} = sortrows(ripple_spikes{day}{epoch}{2,ind});
end

% count only unique neurons firing
ripple_cell_count{day}{epoch} = [cellfun(@unique_neuron_count,ripple_spikes{day}{epoch}(2,:))];
ripple_spikeslist{day}{epoch} = ripple_spikes{day}{epoch};
ripple_barrierlist{day}{epoch}=ripple_barrier;

filename=strcat(f.animal{2},f.animal{1},'ripple_cell_count',dsz,num2str(day));
save(filename,'ripple_cell_count','ripple_spikeslist','ripple_barrierlist');

out.ripple_cell_count = [cellfun(@unique_neuron_count,ripple_spikes{day}{epoch}(2,:))];
out.ripple_spikeslist = ripple_spikes{day}{epoch};
out.ripple_barrierlist=ripple_barrier;

% This will plot the cells of each individual ripple
%cellfun(@(x)(plot_ripple_neurons(x,day,epoch)),ripple_spikes{day}{epoch}(2,:));


disp('done');
%out = 1;

% count the number of unique cells per ripple (limit 1000 cells in
% ripple_spikes)
function num = unique_neuron_count(ripple_spikes)
if(~isempty(ripple_spikes))
    uniqueids = ripple_spikes(:,8).*1000 + ripple_spikes(:,9);
    num = length(unique(uniqueids));
else
    num = 0;
end

% plot each ripple, subplot_count is a global variable that must be set by
% the filter script that calls the entire ripdata filter
function plot_ripple_neurons(ripple_spikes,day,epoch)
global subplot_count;
if unique_neuron_count(ripple_spikes) >= 3
    figure(ceil(subplot_count/12));
    subplot(6,2,mod(subplot_count,12)+1)
    subplot_count = subplot_count + 1;
    uniqueids = ripple_spikes(:,8).*10 + ripple_spikes(:,9);
    plot(uniqueids, 'o');
    title(sprintf('Day %d, Epoch %d',day,epoch));
end