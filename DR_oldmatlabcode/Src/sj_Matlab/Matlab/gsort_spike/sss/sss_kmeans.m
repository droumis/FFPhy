function spikes = sss_kmeans(spikes, options, nch, Subset)

% SS_KMEANS  K-means clustering.
%     SPIKES = SS_KMEANS(SPIKES) takes and returns a spike-sorting object SPIKES.
%     It k-means clusters the (M x N) matrix SPIKES.WAVEFORMS and stores the
%     resulting group assignments in SPIKES.OVERCLUSTER.ASSIGNS, the cluster
%     centroids in SPIKES.OVERCLUSTER.CENTROIDS, and the mean squared error in
%     SPIKES.OVERCLUSTER.MSE.  The W, B, and T matrices (within-cluster, between-
%     cluster, and total) covariance matrices are in SPIKES.OVERCLUSTER.W, etc.
%
%     K-means algorithm is an EM algorithm that finds K cluster centers in
%     N-dimensional space and assigns each of the M data vectors to one of these
%     K points.  The process is iterative: new cluster centers are calculated as
%     the mean of their previously assigned vectors and vectors are then each 
%     reassigned to their closest cluster center.  This finds a local minimum for
%     the mean squared distance from each point to its cluster center; this is
%     also a (local) MLE for the model of K isotropic gaussian clusters.
% 
%     This method is highly outlier sensitive; MLE is not robust to the addition
%     of a few waveforms that are not like the others and will shift the 'true'
%     cluster centers (or add new ones) to account for these points.  See
%     SS_OUTLIERS for one solution.
% 
%     The algorithm used here speeds convergence by first solving for 2 means,
%     then using these means (slightly jittered) as starting points for a 4-means
%     solution.  This continues for log2(K) steps until K-means have been found.
%     Clusters of size one are not allowed and are lumped into the nearest
%     non-singleton cluster, with the result that occasionally fewer than K 
%     clusters will actually be returned.
%
%     SPIKES = SS_KMEANS(SPIKES, OPTIONS) allows specification of clustering
%     parameters.  OPTIONS is a structure with some/all of the following fields
%     defined.  (Any OPTIONS fields left undefined (or all fields if no OPTIONS
%     structure is passed in) uses its default value.)
% 
%         OPTIONS.DIVISIONS (default: round(log2(M/400))) sets the desired number
%               of clusters to 2^DIVISIONS.  The actual number of clusters may be
%               less than this number of singleton clusters are encountered.
%         OPTIONS.REPS (default: 1) specifies the number of runs of the full
%               k-means solution.  The function will return the assignments that
%               resulted in the minimum MSE.
%         OPTIONS.REASSIGN_CONVERGE (default: 0) defines a convergence condition
%               by specifying the max number of vectors allowed to be reassigned
%               in an EM step.  If <= this number of vectors is reassigned, the
%               this condition is met.
%         OPTIONS.MSE_CONVERGE (default: 0) defines a second convergence condition.
%               If the fractional change in mean squared error from one iteration
%               to the next is smaller than this value, this condition is met.
%
%     NOTE: Iteration stops when either of the convergence conditions is met.
%  
% References:
%     Duda RO et al (2001).  _Pattern Classification_, Wiley-Interscience
%
% Last Modified: sbm, 10/03/03

% Undocumented option: OPTIONS.PROGRESS (default: 1) determines whether the progress
%                                 bar is displayed during the clustering.

starttime = clock;

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') || (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
end

if nargin<3,
    nch = 4;
end

if ((nargin < 4) || isempty(Subset)),
    Subset = 1:size(spikes.waveforms,2)/nch;
end

%%%%%%%%%% CONSTANTS
%waves = spikes.waveforms;         % using a reference without structure notation is better for R13 acceleration
Nsamples = size(spikes.waveforms,2)/nch;

if Nsamples > length(Subset)
    waves=[];
    for i=1:nch,
        waves = [waves, spikes.waveforms(:,Nsamples*(i-1) + Subset)];
    end
    spikes.waveforms = waves;
    clear waves;
end

[M,N] = size(spikes.waveforms);
target_clustersize = 500;
jitter = meandist_estim(spikes.waveforms) / 100 / N;        % heuristic

%%%%%%%%%% DEFAULTS

opts.divisions = round(log2(M / target_clustersize));  % power of 2 that gives closest to target_clustersize, but in [4..7]
opts.divisions = max(min(opts.divisions, 8), 4);       %     restrict to 16-256 clusters; heuristic.
opts.memorylimit = M;                   % all spikes at once
opts.memorylimit = 300;                 % 
opts.reps = 1;   
opts.reassign_converge = 0;   % stop only when no point are reassigned ...
opts.reassign_rough = round(0.005*M);   % (no need at all to get full convergence for intermediate runs) 
opts.mse_converge = 0;                  % ... and by default, we don't use the mse convergence criterion
opts.progress = 1;
opts.splitmore = 1;

if (nargin > 1)
	supplied = lower(fieldnames(options));   % which options did the user specify?
	for op = 1:length(supplied)              % copy those over the defaults
        if (version('-release') < 13)        % annoyingly, pre-R13 Matlab doesn't do dynamic field names, ...
            opts = setfield(opts, supplied(op), getfield(options, supplied(op)));  % so we use an older syntax
        else
            opts.(supplied{op}) = options.(supplied{op});  % this is the preferred syntax as of R13 --
        end                                                %   we include it b/c 'setfield' is deprecated
	end
end

if isempty(opts.divisions)
    opts.divisions = round(log2(M / target_clustersize));  % power of 2 that gives closest to target_clustersize, but in [4..7]
    opts.divisions = max(min(opts.divisions, 8), 4);       %     restrict to 16-256 clusters; heuristic.
end
% if ((nargin < 2) || isempty(options.reps))
%     opts.reps = 1;                % just one repetition
% end
% if ((nargin < 2) || isempty(options.reassign_converge))
%     opts.reassign_converge = 0;   % stop only when no point are reassigned ...
% end
% if ((nargin < 2) || isempty(options.mse_converge))
%     opts.mse_converge = 0;        % ... and by default, we don't use the mse convergence criterion
% end
% if ((nargin < 2) || isempty(options.progress))
%     opts.progress = 1;
% end

%%%%%%%%%% CLUSTERING
normsq = sum(spikes.waveforms.^2, 2);
assigns = ones(M, opts.reps);
mse = Inf * ones(1, opts.reps);
% clear fast_kmeans;

for rep = 1:opts.reps                                 % TOTAL # REPETITIONS
	
	centroid = mean(spikes.waveforms, 1);  % always start here	
	for iter = 1:opts.divisions                       % # K-MEANS SPLITS
		itercounter(iter) = 0;
        oldmse = Inf;
		oldassigns = zeros(M, 1);
		assign_converge = 0;
		mse_converge = 0;
        
        if (iter == opts.divisions), progress = opts.progress;  reassign_criterion = opts.reassign_converge; 
		else,                        progress = opts.progress;              reassign_criterion = opts.reassign_rough;
		end
        
		if ((iter==opts.divisions) && opts.splitmore)  % do some extra splitting on clusters that look too big
            toobig = find(clustersizes > 4*target_clustersize);
            splitbig = centroid(toobig,:) + randn(length(toobig),size(centroid,2));
            centroid = [centroid; splitbig];
        end

		centroid = [centroid; centroid] + jitter * randn(2*size(centroid, 1), size(centroid,2)); % split & jitter
        clustersizes = ones(1,size(centroid,1));  % temporary placeholder to indicate no empty clusters at first
        
		while (~(assign_converge | mse_converge))     % convergence?
			% Vectorized distance computation:
			%         dist(x,y)^2 = (x-y)'(x-y) = x'x + y'y - 2 x'y,
			% which can be done simultaneously for all spikes.
			numclusts = size(centroid,1);  % this might not be 2^iter
			centernormsq = sum(centroid.^2, 2);
			clustdistsq = [normsq * ones(1,numclusts)] + [ones(M, 1) * centernormsq']; % faster than repmat (in R13 anyway)
			%clustdistsq = repmat(normsq, [1,numclusts]) + %repmat(centernormsq', [M,1]);
            clustdistsq = clustdistsq - (2 * spikes.waveforms * centroid');  % don't bother taking sqrt; doesn't affect min

			%%%%% E STEP
			% Write cluster assignments by finding where the minimum distance^2 occurs
			[bestdists, assigns(:,rep)] = min(clustdistsq, [], 2);
			
			% Prevent singleton clusters because they screw up later calculations.
			clustersizes = hist(assigns(:,rep), 1:numclusts);
			while (any(clustersizes == 1))
				singleton = find(clustersizes == 1);
				singleton = singleton(1);                          % pick a singleton ...
				singlespike = find(assigns(:,rep) == singleton);   % ... and its lone member ...
				clustdistsq(singlespike, singleton) = Inf;         % ... and find its next best assignment
				[bestdists(singlespike), assigns(singlespike, rep)] = min(clustdistsq(singlespike,:));
				
				clustersizes(singleton) = 0;                       % bookkeeping for cluster counts
				clustersizes(assigns(singlespike, rep)) = clustersizes(assigns(singlespike, rep)) + 1;
			end
			
			%%%%% M STEP
			% Recompute cluster means based on new assignments . . .
			for clust = 1:numclusts
				if (clustersizes(clust) > 0)
					members = find(assigns(:,rep) == clust);
					centroid(clust,:) = mean(spikes.waveforms(members,:), 1);
				end
			end
			centroid(find(clustersizes == 0),:) = [];
			
			%%%%% Compute convergence info
			mse(rep) = mean(bestdists);
			mse_converge = ((1 - (mse(rep)/oldmse)) <= opts.mse_converge);   % fractional change
			oldmse = mse(rep);
			
			changed_assigns = sum(assigns(:,rep) ~= oldassigns);
			assign_converge = (changed_assigns <=  reassign_criterion);   % num waveforms reassigned
            if (progress)
                progressBar(((M - changed_assigns)/M).^10, 5, ['Calculating ' num2str(2^iter) ' means.'] ); % crude ...
            end
			oldassigns = assigns(:,rep);	
            itercounter(iter) = itercounter(iter) + 1;
		end
    end
    
    % Finally, reassign spikes in singleton clusters to next best ...
    assignlist = assigns(:,rep);
	junkclust = find([clustersizes <= 1]);   junkspikes = ismember(assignlist, junkclust);
	centroid(junkclust,:) = Inf;  % take these out of the running & recompute distances
	dist_to_rest = pairdist(spikes.waveforms(junkspikes,:), centroid, 'nosqrt');
	[bestdists(junkspikes), assignlist(junkspikes)] = min(dist_to_rest, [], 2);
	mse(rep) = mean(bestdists);
	
	assigns(:,rep) = assignlist;
end

spikes.overcluster.iteration_count = itercounter;
%spikes.overcluster.randn_state = randnstate;
%spikes.overcluster.rand_state = randnstate;

% Finish up by selecting the lowest mse over repetitions.
[bestmse, choice] = min(mse);
spikes.overcluster.assigns = sortAssignments(assigns(:,choice));
spikes.overcluster.mse = bestmse;
spikes.overcluster.sqerr = bestdists;

% We also save the winning cluster centers as a convenience
numclusts = max(spikes.overcluster.assigns);
spikes.overcluster.centroids = zeros(numclusts, N);
for clust = 1:numclusts
	members = find(spikes.overcluster.assigns == clust);
    spikes.overcluster.centroids(clust,:) = mean(spikes.waveforms(members,:), 1);
end

% And W, B, T matrices -- easy since T & B are fast to compute and T = W + B)
spikes.overcluster.T = cov(spikes.waveforms);             % normalize everything by M-1, not M
spikes.overcluster.B = cov(spikes.overcluster.centroids(spikes.overcluster.assigns, :));
spikes.overcluster.W = spikes.overcluster.T - spikes.overcluster.B;

% Finally, assign colors to the spike clusters here for consistency ...
cmap = jetm(numclusts);
spikes.overcluster.colors = cmap(randperm(numclusts),:);

if (opts.progress), progressBar(1.0, 1, ''); end
spikes.tictoc.kmeans = etime(clock, starttime);