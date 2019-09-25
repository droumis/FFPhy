function spikes = sss_energy(spikes)

% SS_ENERGY  Interface energy based cluster similarity computation.
%     SPIKES = SS_ENERGY(SPIKES) adds an upper triangular interface energy
%     matrix to a spike-sorting object in SPIKES.HIERARCHY.INTERFACE_ENERGY.
% 
%     The energy similarity matrix is calculated by applying an exponential
%     decay to all pairwise euclidean distances between waveforms from two
%     clusters (or within a single cluster for intra-cluster energy) and 
%     summing these distances.
%
%     The calculation ignores the energy due to the zero distance between
%     points and themselves; this removes a dependence of the density on
%     the absolute size of the cluster.  As a result, singleton clusters
%     do not have a well-defined energy and will error out.
%  
%     When each entry is normalized by the number of distinct contributing
%     pairs (Na*Nb for off diagonal entries and Na*(Na-1)/2 on the diagonal),
%     it approximates the fraction of pairs in a given cluster whose distance
%     is not much greater than the length constant of the exponential and thus
%     provides an estimate of local density.  This function does not, however,
%     normalize SPIKES.HIERARCHY.INTERFACE_ENERGY, since the normalized form is
%     inconvenient during cluster aggregation.  The normalization can readily
%     be done, however, with
%          normalize = ((numpts * numpts') - diag(numpts));
%          normalize = normalize - diag(0.5*diag(normalize));
%          normalized_energy = interface_energy ./ normalize;
%     where 'numpts' is a vector of cluster sizes.
%
%     The unnormalized energy matrix can be updated during aggregation without
%     the need to recompute it from scratch.  The intra-cluster energy E(AB,AB)
%     of a cluster AB formed by aggregating clusters A and B is given by
%              E(AB,AB) = E(A,A) + E(B,B) + E(A,B)
%     and the inter-cluster energy between any cluster C and an aggregate AB is
%                 E(AB,C) = E(A,C) + E(B,C)
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%
% Last Modified: sbm, 8/2103

starttime = clock;

%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') | (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
elseif (~isfield(spikes, 'overcluster'))
    error('SS:overcluster_not_computed', 'The data must be overclustered before computing energy');
end
numclusts = length(unique(spikes.overcluster.assigns));
waves = spikes.waveforms;

%%% Take subset of waveforms within index 6:24 from each electrode for energy calculation
%% NEED A SPIKES.NCH FIELD TO DO THIS RIGHT %%
% waves = [waves(:,6:24), waves(:,38:56), waves(:,70:88), waves(:,102:120)]; 


%%%%% PREPARE SOME INFORMATION
normsqr = sum(waves.^2,2);
pts = cell(numclusts,1);    % collect spike indices for each cluster
for clust = 1:numclusts   
    pts{clust} = find(spikes.overcluster.assigns == clust);
end
numpts = cellfun('length', pts);
if (any(numpts < 2))
    error('SS:energy_ill_defined', 'Clusters with fewer than 2 points do not have a defined energy.');
end

%%%%% HEURISTIC DISTANCE SCALE suggested in Fee & seems to work.  The parameter is not too sensitive.
scale = sqrt(sum(diag(spikes.overcluster.W))) ./ 10;

%%%%% PREPARE TO LOOP
total = (numclusts^2 + numclusts) / 2;
k = 1;
progressBar(0, max(floor(total/100),1), 'Computing Interaction Energies . . .')
interface_energy = zeros(numclusts);

%%%%% PAIRWISE DISTANCES LOOP
assigns = spikes.overcluster.assigns;
members = repmat(NaN, [max(numpts), numclusts]);
for clust = 1:numclusts
    local = find(assigns == clust);
    members(1:length(local), clust) = find(assigns == clust);
end
for clust1 = 1:numclusts
    for clust2 = clust1:numclusts   % clust2 starts at clust1 so we get intra- too 
        % Vectorized distance computation inlined for efficiency:
        %         dist(x,y)^2 = (x-y)'(x-y) = x'x + y'y - 2 x'y
        dists = repmat(normsqr(pts{clust1}), [1,numpts(clust2)]);
        dists = dists + repmat(normsqr(pts{clust2})', [numpts(clust1),1]);
        dists = dists - (2 * waves(pts{clust1},:) * waves(pts{clust2},:)');
        dists = sqrt(abs(dists));
        interface_energy(clust1,clust2) = sum(exp(-dists(:)./scale));
 
%         members1 = members(1:numpts(clust1),clust1);
%         members2 = members(1:numpts(clust2),clust2);
%         normsqr1 = normsqr(members1);
%         normsqr2 = normsqr(members2);
%         dists = normsqr1(:, ones(length(members2), 1));
%         dists = dists + normsqr2(:, ones(length(members1), 1))';
%         dists = dists - (2 * waves(members1,:) * waves(members2,:)');
%         dists = sqrt(abs(dists));
%         interface_energy(clust1,clust2) = sum(exp(-dists(:)./scale));
 
        k = k + 1;
        progressBar(k/total);
    end
    clust1
end

%%%%% CORRECTION TERMS
% The energy matrix so far includes a contribution in the intra-cluster
% energies that is not found in the inter-cluster energies; namely, the
% computation of   sum_(all x) sum_(all y) e^(-dist/scale)   for
% intra-cluster energy includes cases where x == y (so dist == 0).
interface_energy = interface_energy - diag(numpts);     % So subtract this out.

% Also, we've double counted pairs in the intra-energy case, since dist(a,b)
% and dist(b,a) are not treated as distinct;
interface_energy = interface_energy - diag(0.5*diag(interface_energy));

%%%%% FINISH UP
spikes.hierarchy.interface_energy = interface_energy;
spikes.tictoc.energy = etime(clock, starttime);
