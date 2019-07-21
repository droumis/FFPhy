
% get performance, change, speed, per ripple
% load lfpstack to get the riptimes and dayeps
% load behavestate
% loop through the dayeps
% filter out any ripple before first and after the last well visit
% get the behavestate statespace allbound, outbound, inbound that are in
% this dayep 

% fuck i still need to exclude rips after the last well visit

load_stack = 0;

get_ripstate = 0;
load_ripstate = 0;
get_designmat = 1;
load_designmat = 0;

load_raw_pwr = 0;
run_regression = 0;
run_multiple_regress = 0;

Fp.animals = {'D10'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz'};

Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
Fp.uselfptype = 'eeg';
Fp.useripstates = {'onlywdays','rewarded', 'unrewarded', 'inbound' , 'outbound'};
pconf = paramconfig;
me = animaldef('Demetris');

%% load lfpstack
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% make, get ripstates.. 

if get_ripstate
    ripstate = getStateFilters(lfpstack);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end

% right now it's a one hot encoded condition filter
% need to also include design matrix of continuous variables
% need to also include a design matrix of discrete/categorical variables
% for now i'll just make the design matrix a seperate thing.. 
%% make design matrix
if get_designmat
    designmat = getDesignMatrix(lfpstack);
    save_data(designmat, [pconf.andef{2}, 'designmat/'], 'designmat_wtrack');
end
if load_designmat
    designmat = load_data([pconf.andef{2}, 'designmat/'], 'designmat_wtrack', Fp.animals);
end

%% load raw power
if load_raw_pwr
    rawpwr = struct;
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        wp = getWaveParams(Fp.waveSet);
        rawpwr = load_data(sprintf('%s/analyticSignal/', me{2}), ...
            sprintf('AS_waveSet-%s_%s_power', wp.waveSet, Fp.uselfptype), animal);
    end
end
%% run single regression
if run_regression
    % load the raw power and do a tiedrank of the power across trials (for each time freq point)
    % reshape the power from nt x time x rip x freq to rip x timefreq for each nt
    %     istate = Fp.useripstates{rs};
    %     istateidx = find(strcmp(istate, ripstate.statesetsfields));
    Pp = load_plotting_params({'defaults', 'riptriglfp_perstatefrex_allntrodes'});
    wp = getWaveParams('4-300Hz');
    nti = 5;
    anidx = 1;
    
    voxel_pval = 0.01;
    mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
    mcc_cluster_pval = 0.05;
    
    % get speed, performance, change, per rip
%     rstartrank = tiedrank(ripstate.ripStartTime');
    rstartrank = tiedrank(ripstate.dayeps(:,1)');
    numrips = length(rstartrank);
    
    include_rips = ripstate.statesets(:,1);
    data = permute(squeeze(rawpwr(anidx).pwr(nti, :, find(include_rips),:)),[2 1 3]);
    data = trim2win(data, Fp.srate, Pp.pwin, 'dsamp', rawpwr(anidx).wp.dsamp);
    nTimepoints = size(data,2);
    data = permute(data, [3 2 1]);
    powerreshaped = reshape(data, wp.numfrex*nTimepoints,numrips)';
    powerrank = tiedrank(powerreshaped); % tiedrank on dim 1, across rips
    
    realcorrs = (rstartrank*rstartrank')\rstartrank*powerrank;
    realcorrs = reshape(realcorrs,wp.numfrex,nTimepoints);
    
    ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(data,2));
    
    sf = subplot(221);
    
    contourf(ptime,wp.frex,realcorrs,40,'linecolor','none')
    title('realcorrs')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    set(gca,'ydir','normal','yscale','log');
    caxis(sf, 'auto')
    ytickskip = 2:4:wp.numfrex;
    set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
%     title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
%         'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
    yl = ylim;
    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
    
    n_permutes = 1000;
    % initialize null hypothesis matrices
    permuted_rvals  = zeros(n_permutes,numfrex,nTimepoints);
    max_pixel_rvals = zeros(n_permutes,2);
    max_clust_info  = zeros(n_permutes,1);
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:n_permutes
        fake_rt_mapping = rstartrank(randperm(numrips));
        
        % compute t-map of null hypothesis
        fakecorrs = (fake_rt_mapping*fake_rt_mapping')\fake_rt_mapping*powerrank;
        
        % reshape to 2D map for cluster-correction
        fakecorrs = reshape(fakecorrs,wp.numfrex,nTimepoints);
        
        % save all permuted values
        permuted_rvals(permi,:,:) = fakecorrs;
        
        % save maximum pixel values
        max_pixel_rvals(permi,:) = [ min(fakecorrs(:)) max(fakecorrs(:)) ];
    end

    % this time, the cluster correction will be done on the permuted data, thus
    % making no assumptions about parameters for p-values
    for permi = 1:n_permutes
        
        % indices of permutations to include in thresholding at this iteration
        perms2use4distribution = true(1,n_permutes);
        perms2use4distribution(permi) = 0;
        
        % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
        fakecorrsz = squeeze((permuted_rvals(permi,:,:)-mean(permuted_rvals(perms2use4distribution,:,:),1)) ./ std(permuted_rvals(perms2use4distribution,:,:),[],1) );
        fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
        
        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(fakecorrsz);
        max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    end
    
    % now compute Z-map
    zmap = (realcorrs-squeeze(mean(permuted_rvals,1)))./squeeze(std(permuted_rvals));
    
    subplot(222)
    contourf(ptime,frex,zmap,40,'linecolor','none')
    title('Unthresholded Z map')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    set(gca,'ydir','normal','yscale','log');
    caxis(sf, 'auto')
    ytickskip = 2:4:wp.numfrex;
    set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
    %     title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
    %         'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
    yl = ylim;
    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
    
    
    % apply pixel-level corrected threshold
    lower_threshold = prctile(max_pixel_rvals(:,1),    mcc_voxel_pval*100/2);
    upper_threshold = prctile(max_pixel_rvals(:,2),100-mcc_voxel_pval*100/2);
    
    zmapthresh = zmap;
    zmapthresh(realcorrs>lower_threshold & realcorrs<upper_threshold)=0;
    subplot(223)
    contourf(ptime,frex,zmapthresh,40,'linecolor','none')
%     axis square
%     set(gca,'clim',[-4 4],'xlim',[-500 1200])
    title('Pixel-corrected Z map')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    set(gca,'ydir','normal','yscale','log');
    caxis(sf, 'auto')
    ytickskip = 2:4:wp.numfrex;
    set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
    %     title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
    %         'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
    yl = ylim;
    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
    
        
    
    % apply cluster-level corrected threshold
    zmapthresh = zmap;
    % uncorrected pixel-level threshold
    zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(zmapthresh);
    clust_info = cellfun(@numel,clustinfo.PixelIdxList);
    clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);
    
    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    
    subplot(224)
    contourf(ptime,frex,zmapthresh,40,'linecolor','none')
%     axis square
%     set(gca,'clim',[-4 4],'xlim',[-500 1200])
    title('Cluster-corrected Z map')
%     xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    set(gca,'ydir','normal','yscale','log');
    caxis(sf, 'auto')
    ytickskip = 2:4:wp.numfrex;
    set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
    %     title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
    %         'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
    yl = ylim;
    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
    
        
end

    %% multiple regession
if run_multiple_regress

    % define covariates (RT and trial number)
    % load the ripstate (design matrix) and zscore the continuous variables
    % X = [ zscore(rts') o1power' ]';
    X = zscore(ripstates);
    
    % permutation test
    % shuffle the design matrix, save the values
    % apply cluster correction
    
    
    %% Figure 34.6
    
    voxel_pval = 0.01;
    mcc_cluster_pval = 0.05;
    
    % note: try to use 1000 or more permutations for real data
    n_permutes = 1000;
    
    realbeta = (X*X')\X*powerrank';
    realbeta = reshape(realbeta,[2 numfrex nTimepoints]);
    
    % initialize null hypothesis matrices
    permuted_bvals = zeros(n_permutes,2,numfrex,nTimepoints);
    max_clust_info = zeros(n_permutes,2);
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:n_permutes
        
        % randomly shuffle trial order
        fakeX = X(:,randperm(EEG.trials));
        
        % compute beta-map of null hypothesis
        fakebeta = (fakeX*fakeX')\fakeX*powerrank';
        
        % reshape to 2D map for cluster-correction
        fakebeta = reshape(fakebeta,[2 numfrex nTimepoints ]);
        
        % save all permuted values
        permuted_bvals(permi,:,:,:) = fakebeta;
    end
    
    % this time, the cluster correction will be done on the permuted data, thus
    % making no assumptions about parameters for p-values
    for permi = 1:n_permutes
        
        for testi=1:2
            % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
            fakecorrsz = squeeze((permuted_bvals(permi,testi,:,:)-mean(permuted_bvals(:,testi,:,:),1)) ./ std(permuted_bvals(:,testi,:,:),[],1) );
            fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
            % get number of elements in largest supra-threshold cluster
            clustinfo = bwconncomp(fakecorrsz);
            max_clust_info(permi,testi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
        end
    end
    
    
    figure
    for testi=1:2
        
        % now compute Z-map
        zmap = (squeeze(realbeta(testi,:,:))-squeeze(mean(permuted_bvals(:,testi,:,:),1))) ./ squeeze(std(permuted_bvals(:,testi,:,:),[],1));
        
        subplot(2,3,1+(testi-1)*3)
        contourf(tftimes,frex,zmap,40,'linecolor','none')
        axis square
        set(gca,'clim',[-3 3],'xlim',[-500 1200])
        title('Unthresholded Z map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
        
        % apply uncorrected threshold
        zmapthresh = zmap;
        zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
        subplot(2,3,2+(testi-1)*3)
        contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
        axis square
        set(gca,'clim',[-3 3],'xlim',[-500 1200])
        title('Uncorrected thresholded Z map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
        
        % apply cluster-level corrected threshold
        zmapthresh = zmap;
        uncorrected pixel-level threshold
        zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
        % find islands and remove those smaller than cluster size threshold
        clustinfo = bwconncomp(zmapthresh);
        clust_info = cellfun(@numel,clustinfo.PixelIdxList);
        clust_threshold = prctile(max_clust_info(:,testi),100-mcc_cluster_pval*100);
        
        % identify clusters to remove
        whichclusters2remove = find(clust_info<clust_threshold);
        
        % remove clusters
        for i=1:length(whichclusters2remove)
            zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
        end
        
        subplot(2,3,3+(testi-1)*3)
        contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
        axis square
        set(gca,'clim',[-3 3],'xlim',[-500 1200])
        title('Cluster-corrected Z map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    end
end