function [spikes,draw] = sss_spksort(spkfile, parameters)


%% Demonstration code for the spike sorter derived from the Fee et al.
%% (J. Neurosci Meth (1996) 69: 175-88) algorithm.  Requires Matlab
%% version 6.1 or later and the Statistics Toolbox.
%%                                          Last Modified: 8/22/04
%%                                          samar mehta, sbmehta@ucsd.edu
%%


% NOTE: PARAMETERS CAN BE MUCH BETTER HANDLED BY A VARARGIN OPTION, BUT I AM
% DOING IT RIGHT NOW BY CHECKING THAT EACH PARAMETER HAS A VALUE BEFORE CALLING
% THE APPROPRITATE FUNCTION:
%                            Shantanu Jadhav, shantanu@ucsd.edu


if (nargin < 2)
   parameters.outliers_a =[];
   parameters.kmeans_options.divisions=[];
   parameters.reint_out=[]; parameter.tmin=[]; parameters.tref=[]; parameters.cutoff = [];
end


% LOAD OUTPUT FILE FROM AUTOSPIKECUT/MAKEONEFILE & THEN RUN

load(spkfile);
draw.raw=spikes.waveforms;
   
%%%%% ALIGNING SPIKES
%% Noise on the electrode can jitter the exact time at which threshold crossing
%% occurs.  If you are using this method to extract your spikes, this can be
%% a significant source of variability for spikes from the same neuron.  The
%% next step 'de-jitters' the data by aligning all spikes on the same sample.
%% (Note: this requires THRESHT and THRESHV to be defined as described above).
%% The data are then replotted just as they were in figure 1.  Put the plots
%% side by side.  The density plot for the dejittered data should look more tight.

spikes = sss_dejitter(spikes);
draw.dejitter = spikes.waveforms;
%     h=figure; colormap hot;
%     subplot(2,1,1); plot(spikes.waveforms'); axis tight; title('Centered Data');
%     subplot(2,1,2); hist2d(spikes.waveforms); clim = get(gca, 'CLim'); %set(h, 'CLim', clim);

%%%%% REMOVING OUTLIERS
%% The density estimation techniques used in spike sorting routines are typically
%% not robust -- that is, they are sensitive to outliers -- so we need to remove
%% outliers.  'outliers' are spikes that look so much unlike the other events in
%% the data that they cannot be sorted reliably.  This does not mean that they are
%% not interesting (they often include overlapping spikes/doublets in addition to
%% electrical artifact), it just means that you'll have to look at them by hand.
%% The following code removes the worst offenders.  See the help for 'ss_outliers'
%% if you want to be more/less conservative.  (We don't replot density because
%% this won't change much since outliers are only a small percent of the data).

spikes = sss_outliers(spikes, parameters.outliers_a);

draw.nooutliers = spikes.waveforms;
draw.outliers = spikes.outliers.waveforms;

    % figure(5); 
    % plot(spikes.waveforms'); axis tight; title('Centered Data w/ Outliers Removed');


%%%%% INITIAL CLUSTERING
%% The Fee algorithm deals with possibly non-Gaussian data (e.g., bursting neurons)
%% by doing the sorting in two steps.  The first step fits many local Gaussian
%% spheres to the data to identify groups of spikes with similar shapes; these will
%% later be combined into spike assignments.  This two step procedure is a good place
%% to do a sanity check; do the results of the local clustering look like the
%% algorithm is capturing local density?  The following code plots the waveforms
%% (now colored according to local similarity) and the type of height-width that
%% is often used for manual sorting (colored similarly).

    spikes = sss_kmeans(spikes, parameters.kmeans_options);

    % figure(6);  set(gcf, 'Renderer', 'OpenGL');
    % clustshow(spikes, spikes.overcluster.assigns, [1:32]); subplot(2,1,1); title('Local Clusters');



%%%%% FINAL CLUSTERING
%% The local clusters are now combined into larger groups corresponding (with luck)
%% to neurons.  This step is typically the most time consuming (especially for
%% larger data sets).  After the aggregation, the final results are summarized in
%% a series of plots.
%% Figure 5 contains one row of graphs for each final cluster.  The left plot is
%% a density image and the two right plots are interspike interval histograms
%% (at two time scales).  If the clustering is good, there should be few events
%% with interspike intervals less than 2 msec; the ISI score displayed in the
%% middle column reflects this (smaller numbers are better).  
%% Figure 6 looks like Figure 4 from above, but recolored after aggregation.
%% The bottom plot also contains a legend matching the colors to the cluster #
%% labels in Figure 5.
%% Figure 7 shows how the local clusters were combined into the final clusters.

pack;

spikes = sss_energy(spikes);

%save tempenergy
file1 = [spkfile '_energy'];
save(file1);

agg.reint_out=parameters.reint_out; agg.tmin = parameters.tmin; agg.tref = parameters.tref; agg.cutoff=parameters.cutoff
spikes = sss_aggregate(spikes, agg);

figure(99); colormap hot;
showclust(spikes, spikes.hierarchy.assigns);
draw.sorted = spikes.waveforms;
    
newassigns = spikes.hierarchy.assigns;

    % figure(10); ssg_databrowse2d(spikes);
    % figure(11); ssg_databrowse3d(spikes);

% FIGURES LIKE OLD METHOD

newplots=newassigns;

x_isi=20; tmin=0.001;
arfp = tmin*1000;
skip = [1:4:41]; skipp= skip+1;
clr = ['y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r'];

stimonset=200; window=50; 
sweeps = (spikes.fstimes-mod (spikes.fstimes, 10000))/10000;

sss_rawrastersisi (spikes.waveforms, spikes.fstimes, spikes.hierarchy.assigns, clr, arfp, skip, skipp, spikes.stimonset, x_isi, 7, spikes.sweepd, spikes.nsweeps);
%sss_rawwaves (spikes.waveforms, spikes.fstimes, spikes.hierarchy.assigns, clr, arfp, skip, 8);

%--------------------------------------------------------------------------

%   DO NOT REALLY NEED THIS: INFACT, ADD THIS AT THE END OF SAVE COMMAND IN THE GUI: CAN ALSO THEN ADD COMMENT THERE ABOUT S,M,U  

% lb = unique (newassigns);
% for lp_lb = 1: length (lb)
%     % First row = Cluster labels
%     spikes.clusters {1,lp_lb} = (['Cluster' num2str(lb(lp_lb))]);
%     % Second row = waveform voltage values
%     spikes.clusters {2,lp_lb} = fwaves ((newassigns==lb(lp_lb)),:);
%     % Third row = spiketimes
%     spikes.clusters {3,lp_lb} = fstimes ((newassigns==lb(lp_lb)),:);
%     % Fourth row = comment of cluster being single, multi unit (or unclear)
% end

%------------------------------------------------------------------------


%% (note: For neurons with relatively high firing rates, the algorithm can watch
%%         for refractory periods to determine when to stop the aggregation.  Since
%%         this does not work well for neurons with very low firing rates, it can
%%         sometimes make mistakes towards the end of the aggregation, joining
%%         clusters that should not be joined or vice versa.  A person looking at
%%         figures 5-7 can usually see this fairly quickly.  The problem can be
%%         corrected using the following types of commands:)
%% spikes = merge_clusters(spikes, clusternumber1, clusternumber2);
%% spikes = split_cluster(spikes, clusternumber);
%%         After you do this, replot figures 5-7 above and see if things look better).

%%%%% Getting the cluster assignments
%% The vector of cluster assignments can be found in the following location:
%%                          spikes.hierarchy.assigns

%save tempsort;
file = [spkfile '_sort'];
save(file, 'spikes');
%save tempsort;