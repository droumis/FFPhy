function spikes = sss_aggregate(spikes, parameters)
% SS_AGGREGATE  ISI-restricted heirarchical cluster aggregation.
%     SPIKES = SS_AGGREGATE(SPIKES) takes and returns a spike-sorting object
%     SPIKES after aggregating clusters (requires a previous overclustering
%     and an interface energy matrix calculation).  The aggregation tree is
%     stored in SPIKES.HIERARCHY.TREE and the new assignments are stored in
%     SPIKES.HIERARCHY.ASSIGNS.
%
%     The algorithm computes a similarity/connection matrix using the interface
%     energy.  It then chooses the cluster pair with the highest connection
%     strength, aggregates them, recalculates connection strengths, and then
%     repeats the process.  Cluster aggregation is contingent on passing an
%     interspike interval (ISI) test if SPIKES.SPIKETIMES is defined; if this
%     test is not passed, the pair is not aggregated and aggregation continues.
%     Aggregation stops when either (1) all remaining pairs fail the ISI test
%     or (2) the connection strength drops below a (heuristic) cutoff of 0.01.
%
%     The SPIKES.HIERARCHY.TREE output is a matrix describing the aggregation.
%     Each aggregation step entry produces a row, listed in the order they
%     were performed.  The first two columns are the indices of the clusters
%     that were aggregated; the index assigned to the aggregate for future
%     entries is the lower of the two original indices.  The third column is
%     the connection strength between the clusters before aggregation and the
%     fourth column is the isi statistic for the aggregate (0 if isi statistics
%     are not being used).
% 
%     After aggregation, outliers that were previously removed are typically 
%     reinserted into the spikes list so that the list aligns with the original
%     (pre-sorted) list.  The outlier waveforms and spike times are thus by
%     default added back to the SPIKES.WAVEFORMS and SPIKES.SPIKETIMES fields
%     respectively, after all other aggregation is complete.  These waveforms
%     are assigned the label 0 in SPIKES.HIERARCHY.ASSIGNS and the other, 
%     non-outlier, spikes are renumbered accordingly.  To prevent this from
%     occuring, pass in 0 as the second argument to this function, i.e.,
%     SPIKES = SS_AGGREGATE(SPIKES, 0);  
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%
% Last Modified: sbm, 10/04/03

starttime = clock;

%cutoff = 0.01;   % arbitrarily stop aggregation when overlap density is < 1% of main cluster density



%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'hierarchy') | ~isfield(spikes.hierarchy, 'interface_energy'))
    error('SS:energy_not_computed', 'An energy matrix must be computed before aggregation.');
elseif (~isfield(spikes, 'overcluster') | ~isfield(spikes.overcluster, 'assigns'))
    error('SS:overcluster_not_computed', 'The data must be overclustered before aggregation.');
elseif (isfield(spikes, 'spiketimes') & ~isfield(spikes, 'Fs'))
    error('SS:bad_stop_condition', 'A sampling frequency Fs must be supplied for an ISI stop condition.');
end
if (isfield(spikes.hierarchy, 'assigns') && any(spikes.hierarchy.assigns == 0))
    error('SS:aggregate_after_outliers', 'Aggregation can not be performed after outliers are reintegrated into the data.');
end
if (nargin < 2)
    reintegrate_outliers = 1;
end

%%%%% PARAMETER CHECKING
if ((nargin < 2) || isempty(parameters.reint_out))
    parameters.reint_out = 1;
    reintegrate_outliers = 1;
end
if ((nargin < 2) || isempty(parameters.tmin))
    parameters.tmin = .001;
    % parameters.tmin = size(spikes.waveforms,2)./spikes.Fs;  % minimum possible time btw spikes
end
if ((nargin < 2) || isempty(parameters.tref))
    parameters.tref = 0.002;
    % parameters.tref = max(0.002, tmin*1.5);                 % crude guess as to refractory period
end
if ((nargin < 2) || isempty(parameters.cutoff))
    parameters.cutoff = 0.01;
    % parameters.tref = max(0.002, tmin*1.5);                 % crude guess as to refractory period
end

tmin=parameters.tmin;
tref=parameters.tref;    
reintegrate_outliers = parameters.reint_out;
cutoff=parameters.cutoff;

%%%%% INITIALIZE A FEW THINGS
assignments = spikes.overcluster.assigns;
interface_energy = spikes.hierarchy.interface_energy;
numclusts = max(assignments);
numpts = full(sparse(assignments, 1, 1, numclusts, 1));
tree = [];
untested = ones(numclusts);    % they're all initially untested

handle_fig = figure;
handle_img = imagesc(untested); colormap gray;

%%%%% AGGREGATE HIGHEST CONNECTION STRENGTHS UNTIL ALL TRIED
while (any(any(triu(untested,1))))   % only the upper triangular half is meaningful
    % compute connection strengths from interface energies
    %   first, normalize energy:
    normalize = ((numpts * numpts') - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
    normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
    norm_energy = interface_energy ./ normalize;
    %   then, compute connection strength
    self = repmat(diag(norm_energy), [1,numclusts]);
    connect_strength = 2 .* norm_energy ./ (self + self');
    connect_strength = connect_strength .* (1-eye(numclusts));  % diag entries <- 0, so we won't agg clusters with themselves

    % Find best remaining pair
    remaining = (untested .* connect_strength);
    best = max(remaining(:));           % highest untested connection strength
    
    if (best < cutoff)   % No point continuing if connection strengths have gotten really lousy
        break;
    end
    
    [clust1 clust2] = find(connect_strength == best);  % who're the lucky winners?
    first = min(clust1(1),clust2(1));   % if we get 2 best pairs, just take the 1st
    second = max(clust1(1),clust2(1)); 
    untested(first,second) = 0;         % mark that we've tried this (in the upper half of 'untested')
    set(handle_img, 'CData', untested); title(['Trying ' num2str(first) ' and ' num2str(second)]); drawnow;
    
    % Is this aggregation allowed?
    if (isfield(spikes, 'fstimes'))  % if we were given spike times, use them ...
        t1 = spikes.fstimes(find(assignments == first))/1000;
        t2 = spikes.fstimes(find(assignments == second))/1000;
		
        [allow, scores] = isiQuality(t1, t2, tmin, 0.010, tref, spikes.Fs);
    else  % ... otherwise, there are no restrictions on aggregation
        allow = 1;
        scores = [0 0 0];
    end

    if (allow)      % Bookkeeping ...        
        % Aggregation subsumes the higher index cluster into the lower.  Start by adding
        % (denormalized) interaction energies for the second (higher index) cluster
        % to those of the first and zeroing the old entries of the second.  Because of the
        % triangular structure of the connection matrix, repeat for both rows and columns ...
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

        % update counts vector
        numpts(first) = numpts(first) + numpts(second);
        numpts(second) = 2;   % leaving this as 2 prevents div by zero during normalization above
        
        % make a tree entry for the aggregation we just performed
        tree = cat(1, tree, [first, second, best, scores(3)]);

        % Now actually change the numbers
        assignments(find(assignments == second)) = first;
        
        % Finally, indicate that potential aggregations between the new cluster and 
        % other (nonempty) clusters are untested while pairs involving clusters that
        % have already been emptied should not be tested.
        untested(first,:) = 1;
        untested(:,first) = 1;
        untested(tree(:,2),:) = 0;
        untested(:,tree(:,2)) = 0;
    end
end
close(handle_fig);

spikes.hierarchy.tree = tree;
spikes.hierarchy.assigns = assignments;

spikes.hierarchy.interface_energyn=interface_energy; %SHANTANU: Keep track of interface energies based on current trees
%spikes.newassigns = spikes.hierarchy.assigns;

if (reintegrate_outliers && isfield(spikes, 'outliers') && ~isempty(spikes.outliers.badinds))
    % First, we make room by putting all of the non-outliers back into their original places
    spikes.waveforms(spikes.outliers.goodinds,:) = spikes.waveforms;
    spikes.spiketimes(spikes.outliers.goodinds,:) = spikes.spiketimes;
    
    spikes.swtimes(spikes.outliers.goodinds,:) = spikes.swtimes;
    spikes.ftimes(spikes.outliers.goodinds,:) = spikes.ftimes;
    spikes.fstimes(spikes.outliers.goodinds,:) = spikes.fstimes;
    
    spikes.hierarchy.assigns(spikes.outliers.goodinds) = spikes.hierarchy.assigns;
    %spikes.newassigns(spikes.outliers.goodinds) = spikes.newassigns;
    
    % Then we fill in the outliers ...
    spikes.waveforms(spikes.outliers.badinds,:) = spikes.outliers.waveforms;
    spikes.spiketimes(spikes.outliers.badinds,:) = spikes.outliers.spiketimes;
    
    spikes.swtimes(spikes.outliers.badinds,:) = spikes.outliers.swtimes;
    spikes.fstimes(spikes.outliers.badinds,:) = spikes.outliers.fstimes;
    spikes.ftimes(spikes.outliers.badinds,:) = spikes.outliers.ftimes;
    
    spikes.hierarchy.assigns(spikes.outliers.badinds) = 0;  % ... and add the '0' label.
    %spikes.newassigns(spikes.outliers.badinds) = 0;
    
    % We'll also want to add the assignments to the 'overcluster' list (this is
    % important for post-clustering splitting).
    spikes.overcluster.assigns(spikes.outliers.goodinds) = spikes.overcluster.assigns;
    spikes.overcluster.assigns(spikes.outliers.badinds) = 0;
    
    spikes = rmfield(spikes, 'outliers');  % don't need this any more -- its redundant.
end

spikes.tictoc.aggregate = etime(clock, starttime);