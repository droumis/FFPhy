function [assignments,spikes,draw] = ssm_tetsort(spkfile, parameters, nch, Spiketype, elec, Subset, saveopt, dodraw)
% 
% [assignments,spikes] = ssm_tetsort('07_030510_tet1_gspike', [], 4, 1, '1', [], 1, 0);
%% Run all steps of gsort algo on tetrode data
% spiketype = 0; % dejiter for negative-going spikes - eg. gspikes
% spiketype = 1; % dejitter for positive going spikes - eg. matspikes
% Whether to save full structure or not - assignments always saved by default in a .clu file
% draw structure if for the gsort GUI

% NOTE: PARAMETERS CAN BE MUCH BETTER HANDLED BY A VARARGIN OPTION, BUT I AM
% DOING IT RIGHT NOW BY CHECKING THAT EACH PARAMETER HAS A VALUE BEFORE CALLING
% THE APPROPRITATE FUNCTION:
%                            Shantanu Jadhav, shantanu@phy.ucsf.edu

if ((nargin < 2) || isempty(parameters))
   parameters.outliers_a =1e-7;
   parameters.kmeans_options.divisions=6; % 2^6=64;
   parameters.reint_out=1; parameters.tmin=0.002; parameters.tref=0.003; parameters.cutoff = 0.01;
end

if (nargin < 3), nch=4; end

if (nargin < 4), Spiketype = 1; end

if (nargin < 5), elec='1'; end

if ((nargin < 6) || isempty(Subset)), 
    Subset = []; 
end

if (nargin < 7), saveopt = 1; end

if (nargin < 8), dodraw = 0; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clufile = [spkfile '.clu.' elec];

load(spkfile);

flag=0;
if (~isfield(spikes, 'waveforms')),
    spikes.waveforms=[]; flag=1;
    for i=1:nch,
        cmd = sprintf('spikes.waveforms = [spikes.waveforms, spikes.waveforms_ch%d];',i); eval(cmd);
    end
else
    spikes.waveforms = [spikes.waveforms];
end


if dodraw == 1, draw.raw=spikes.waveforms; end

%%%%% ALIGNING SPIKES
%% Noise on the electrode can jitter the exact time at which threshold crossing
%% occurs.  If you are using this method to extract your spikes, this can be
%% a significant source of variability for spikes from the same neuron.  The
%% next step 'de-jitters' the data by aligning all spikes on the same sample.
%% (Note: this requires THRESHT and THRESHV to be defined as described above).
%% The data are then replotted just as they were in figure 1.  Put the plots
%% side by side.  The density plot for the dejittered data should look more tight.

maxshift=3;
disp('De-jittering ...');
spikes = ssm_dejitter_fortet2(spikes , maxshift, nch);

if dodraw == 1, draw.dejitter = spikes.waveforms; end

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


% disp('Removing outliers ...');
% spikes = ssm_outliers(spikes, parameters.outliers_a, nch, Subset);
% 
% if dodraw == 1,
%     draw.nooutliers = spikes.waveforms;
%     draw.outliers = spikes.outliers.waveforms;
% end


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

disp('Doing k-means ...');
spikes = sss_kmeans(spikes, parameters.kmeans_options, nch, Subset);
% NOTE - 04/23/10, CHANGING FROM SSM_KMEANS HERE BECAUSE DOESNT WORK ON 64-BIT LINUX MACHINE

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

%pack;

disp('Calculating energy - this can take some time - go get coffee  ...');
spikes = ssm_energy_subset(spikes, nch, Subset);

%save tempenergy
file1 = [spkfile '_autosort_energy'];
save(file1,'-v7.3');

agg.reint_out=parameters.reint_out; agg.tmin = parameters.tmin; agg.tref = parameters.tref; agg.cutoff=parameters.cutoff
disp('Aggregating clusters ...');
spikes = ssm_aggregate(spikes, agg);

%figure(99); colormap hot;
%showclustn(spikes, spikes.hierarchy.assigns);
if dodraw == 1,
    draw.sorted = spikes.waveforms;
end

newassigns = spikes.hierarchy.assigns;
unqassigns = unique(spikes.hierarchy.assigns);


%%% SAVE CLU FILE
assignments = spikes.hierarchy.assigns;
nclu = length(unique(assignments)); % Including 0
fid = fopen(clufile, 'wt');
fprintf(fid, '%3.0f\n', nclu);
fprintf(fid, '%3.0f\n', assignments);
fclose(fid);

%%% Remove spikes.waveforms,
if flag==1,
    spikes = rmfield(spikes,'waveforms');
end

    % figure(10); ssg_databrowse2d(spikes);
    % figure(11); ssg_databrowse3d(spikes);


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

disp('Sorted and saved .clu file! If you chose saveopt=1, I am saving "spikes" structure file as well ...');
if saveopt==1
    file = [spkfile '_auto_gsort'];
end
save(file, 'spikes','-v7.3');
