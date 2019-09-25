function spikes = sss_outliers(spikes, a, reps)

% SS_OUTLIERS  K-means based outlier detection.
%     SPIKES = SS_OUTLIERS(SPIKES) takes and returns a spike-sorting object
%     SPIKES.  It identifies likely outliers in the SPIKES.WAVEFORMS and moves
%     them from the (initially M x N) SPIKES.WAVEFORMS matrix to a (P x N)
%     SPIKES.OUTLIERS.WAVEFORMS (usually P << M), removing them the original
%     SPIKES.WAVEFORMS matrix.  If the SPIKES.OUTLIERS.WAVEFORMS matrix already
%     exists, new outliers are added to the existing matrix.
%
%     Outlier spike timestamps are recorded in SPIKES.OUTLIERS.SPIKETIMES
%     (these event times are removed from SPIKES.SPIKETIMES).  The cell matrix
%     SPIKES.OUTLIERS.WHY stores a string describing the reasons that the
%     corresponding waveforms was marked as an outlier.
% 
%     The SS_OUTLIERS function identifies outliers using an ad hoc heuristic
%     based on a k-means clustering; waveforms that end up very far from their 
%     assigned centroid are likely outliers.   The scale for this determination
%     comes from the average of the covariance matrices of each cluster (i.e.,
%     1/(N-1) times the (N x N) within-group sum of squares matrix).  We take
%     this as an approximation to the noise covariance in outlier-free clusters 
%     and note that if the noise were locally Gaussian, then the waveform Mahalanobis
%     distances to assigned cluster means should be roughly Chi^2 distributed
%     (actually, F-distributed but we ignore the refinement for now).
%     NOTE: Outliers damage the k-means solution; the clustering is not robust to
%            gross violations of its (local) Gaussian assumption.  Repeat clustering  
%            on the cleaned data will thus yield a new solution ... which might
%            uncover further outliers.  This function attempts 3 cluster/clean
%            iterations (rule of thumb; tradeoff btw cleaning and run-time),
%            although it stops after any iteration that does not find an outlier.
%
%     SPIKES = SS_OUTLIERS(SPIKES, CUTOFF) specifies the Chi^2 CDF cutoff.
%     The default value is (1 - 1/M), i.e., spikes with Mahalanobis distance
%     such that their Chi^2 are likely to appear less than once by chance are
%     treated as outliers.
%
%     SPIKES = SS_OUTLIERS(SPIKES, CUTOFF, REPS) performs REPS cluster/clean
%     iterations (default: 3).  To specify REPS while using the default CUTOFF,
%     pass in [] for the second argument.
%
%     The reason given for these outliers in SPIKES.OUTLIERS.WHY is "Poor K-Means".
%
%     NOTE: This function performs a k-means clustering and thus overwrites
%           existing k-means assignments in the SPIKES object.
%
% Last Modified: sbm, 10/03/03

starttime = clock;


%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') | (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
end
times_defined = isfield(spikes, 'spiketimes');
sweeps_defined = isfield(spikes, 'sweep');
[M,N] = size(spikes.waveforms);
if ((nargin < 2) || isempty(a))
    a=0.0001; 
end
cutoff = (1 - (a/M));
if (nargin < 3)
    reps = 3;
end

% Initialize the outliers sub-structure
if (~isfield(spikes, 'outliers'))
    spikes.outliers.waveforms = [];
    spikes.outliers.swtimes = [];
    spikes.outliers.fstimes = [];
    spikes.outliers.ftimes = [];
    
    spikes.outliers.why = [];
    if (times_defined)
        spikes.outliers.spiketimes = [];
        spikes.outliers.swtimes = [];
        spikes.outliers.ftimes = [];
        spikes.outliers.fstimes = [];
    end
    if (sweeps_defined)
        spikes.sweep = [];
        spikes.trial = [];
        spikes.stimulus = [];
        spikes.igorstim = [];
    end
    
    spikes.outliers.goodinds = [1:M]';  % We need these to re-insert the outlier 'cluster'
    spikes.outliers.badinds = [];      % into the waveforms matrix after sorting.
end

for cleaning = 1:reps   % Cluster/clean then rinse/repeat.  3 reps does a good job w/o taking too much time

    progressBar((cleaning-1)./reps, 1, ['Removing Outliers: Pass ' num2str(cleaning)]);
    
    %%%%% PERFORM A K-MEANS CLUSTERING OF THE DATA
    opts.mse_converge = 0.001;        % rough clustering is fine
    opts.progress = 0;
 %   opts.divisions = []; options.reassign_converge=[];options.reps =[];
    spikes = sss_kmeans(spikes, opts);
    
    %%%%% MAHALANOBIS DISTANCES:   (x - x_mean)' * S^-1 * (x - x_mean)
    vectors_to_centers = spikes.waveforms - spikes.overcluster.centroids(spikes.overcluster.assigns,:);
    mahaldists = sum([vectors_to_centers .* (spikes.overcluster.W \ vectors_to_centers')'], 2)';
    
    %%%%% SPLIT OFF OUTLIERS
    bad = find(mahaldists > chi2inv(cutoff, N));   % find putative outliers
    
    if (length(bad) == 0)
        break;              % didn't find anything ... no sense continuing
    else
        %%%%% ADD OUTLIERS TO SS OBJECT AND REMOVE FROM MAIN LIST
        spikes.outliers.waveforms = cat(1, spikes.outliers.waveforms, spikes.waveforms(bad,:));
 %       spikes.outliers.why = cat(1, spikes.outliers.why, repmat({'Poor K-Means'}, [length(bad), 1]));
        spikes.outliers.badinds = cat(1, spikes.outliers.badinds, spikes.outliers.goodinds(bad(:)));
        if (times_defined)
            spikes.outliers.spiketimes = cat(1, spikes.outliers.spiketimes, spikes.spiketimes(bad,:));
            spikes.outliers.swtimes = cat(1, spikes.outliers.swtimes, spikes.swtimes(bad,:));
            spikes.outliers.fstimes = cat(1, spikes.outliers.fstimes, spikes.fstimes(bad,:));
            spikes.outliers.ftimes = cat(1, spikes.outliers.ftimes, spikes.ftimes(bad,:));
            
            spikes.spiketimes(bad,:) = [];
            spikes.swtimes(bad,:) = [];
            spikes.fstimes(bad,:) = [];
            spikes.ftimes(bad,:) = [];
        end
        
        if (sweeps_defined)
            spikes.outliers.sweep = cat(1, spikes.outliers.sweep, spikes.sweep(bad,:));
            spikes.outliers.trial = cat(1, spikes.outliers.trial, spikes.trial(bad,:));
            spikes.outliers.stimulus = cat(1, spikes.outliers.stimulus, spikes.stimulus(bad,:));
            spikes.outliers.igorstim = cat(1, spikes.outliers.igorstim, spikes.igorstim(bad,:));
            
            spikes.sweep(bad,:) = [];
            spikes.trial(bad,:) = [];
            spikes.stimulus(bad,:) = [];
            spikes.igorstim(bad,:) = [];
        end
        
        spikes.outliers.goodinds(bad) = [];
        spikes.waveforms(bad,:) = [];
        spikes = rmfield(spikes, 'overcluster');  % this isn't supposed to be public-visible
    end
end

progressBar(1.0, 1, '');

%%%%% TIMING INFORMATION
spikes.tictoc = rmfield(spikes.tictoc, 'kmeans');
spikes.tictoc.outliers = etime(clock, starttime);
