function [spikes, connect_strength]= sss_update_intenergy_onmerge(spikes, to, from)
%%% Make new interface energy field based on current tree, Also return New Conn Str%%
%%Shantanu: Sep2007%%%%%%%

tree=spikes.hierarchy.tree;
interface_energy=spikes.hierarchy.interface_energy;

first=to; second=from;
interface_energy(first,:) = interface_energy(first,:) + interface_energy(second,:);
interface_energy(second,:) = 0;
interface_energy(:,first) = interface_energy(:,first) + interface_energy(:,second);
interface_energy(:,second) = 0;
interface_energy(second,second) = 1;  % keep self-energy at 1 (we may divide by it later)
% since we added rows & columns, some energy values will have spilled over into the
% lower half of the energy matrix (which must be upper triangular).  The next 2 steps
% recover those values.
overflow = tril(interface_energy, -1);   % spillover below diagonal
interface_energy = interface_energy + overflow' - overflow;  % reflect above diagonal

spikes.hierarchy.interface_energyn=interface_energy;

%%%%%%%% Calculations for normalization %%%%%

rawassigns= spikes.overcluster.assigns; rawassigns(find(rawassigns==0))=[];
rawclusts = unique(rawassigns); 
currassigns= spikes.hierarchy.assigns; currassigns(find(currassigns==0))=[];
currclusts = unique(currassigns); 
numclusts =length(rawclusts);

% currassigns=spikes.hierarchy.assigns; currassigns(find(currassigns==0))=[];
% currclusts = unique(currassigns); 

%% UPDATE NPTS FOR CLUSTERS
for i=1:length(rawclusts), 
    numpts(i) = length(find(currassigns==rawclusts(i))); %% If ant rawclust has been merged, it auto goes to zero size in currassigns
end
numpts(find(numpts==0))=2; %% For clusters which have been merged, avoid getting zero size
numpts=numpts';

normalize = ((numpts * numpts') - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
norm_energies = interface_energy ./ normalize;

self = repmat(diag(norm_energies), [1,numclusts]);
connect_strength = 2 .* norm_energies ./ (self + self');
connect_strength = connect_strength .* (1-eye(numclusts));  % diag entries <- 0, so we won't agg clusters with themselves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%