
% get performance, change, speed, per ripple
% fuck i still need to exclude rips after the last well visit

% now that i have the design matrix of continuous vars per rip.. run a
% timefreq corr with the design vars.. use multuple regression?
% plot a full all nts time freq plot of each varcorr zmap per an..
% i want the zmap overlayed with the pixel,cluster corrected sig testing contours



load_stack = 0;

get_ripstate = 0;
load_ripstate = 0;
get_designmat = 0;
load_designmat = 0;

load_raw_pwr = 0;
run_regression = 1;
save_regression = 1;
load_regression = 0;
plot_regression = 0;
run_TFVarcorr = 0;
plot_corrheatmap = 0;
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
%% run single regression % make this into a new function
% the function should return the corrzmapmasked in dims: freq, time,
% ntrode, expvar.. so 25 x 1501 x 30 x 10 cube

% then 1 plot per expvar (timefreqmap/nt) for each animal.. so 10 plots per
% animal

% then once i have 10 plots for each animal.. i need print them.. identify the
% nt's ... and the var/s ***********
% for continuous vars in the design matrix.. the corr value that is the
% average within a large-small timefreq bin per animal.. is the timefreqvar
% corr for that animal.. then how to combine across animals?? 


% Group level analysis for discrete conditions.. day, epoch, outbound, inbound,
% rewarded, unrewarded:
% use 2b strategy from chapter 35.. for one animal, be
% exploratory to find the strongest relationships (timefreq, ntarea x expvar) 
% then use the rest of the animals to test specific. use large timefreq bin
% and take smaller bin within it centered on the timefreq peak within the
% large timefreq bin.. then get the condition average within that small
% window as the value for that animal,nt,condition,timefreqtype

% so it's enough to get an 2squared and p value per animal (per
% var/ntareatfbox).. so once i pick the nt and tfbox largesmall bounds.. i
% can extract the mean within the small tfbox across rips and get an r2 and
% p for the data vec vs the var vec.. 

%take a design matrix and data matrix 
%% TF x expVars Regression and Zmap
if run_regression
    % load the raw power and do a tiedrank of the power across trials (for each time freq point)
    % reshape the power from nt x time x rip x freq to rip x timefreq for each nt
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        PV(ian).animal = animal;
        PV(ian).designMat = designmat(ian);
        PV(ian).ntrodes = ntrodes;
        PV(ian).expvars = designmat(ian).expvars;
        PV(ian).frequency = rawpwr(ian).frequency;
        for iv = 1:length(designmat(ian).expvars)
            fprintf('expvar %d\n', iv);
            expvarrank = designmat(ian).dm(:,iv)'; 
            numrips = length(find(~isnan(expvarrank)));
            expvarrank = expvarrank(~isnan(expvarrank));
%             include_rips = ones(length(expvarrank),1);
            include_rips = find(~isnan(expvarrank));
            for nti = 1:length(ntrodes)
                fprintf('nt %d \n',nti);

        % get in,iv tfcorr for design vs data
%         [PV(ian).data(in,:,:,iv), PV(ian).zmap(in,:,:,iv), PV(ian).thresh(in,:,:,:)] = ...
%             regress(designmat(ian).data(:,iv), rawpwr(ian).data(nt,win1:win2,:));

                data = permute(squeeze(rawpwr(ian).pwr(nti, :, include_rips,:)),[2 1 3]);
                data = trim2win(data, Fp.srate, Pp.pwin, 'dsamp', rawpwr(ian).wp.dsamp);
                nTimepoints = size(data,2);
                data = permute(data, [3 2 1]);
                powerreshaped = reshape(data, wp.numfrex*nTimepoints,numrips)';
                powerrank = tiedrank(powerreshaped); % tiedrank on dim 1, across rips
                tmp = (expvarrank*expvarrank')\expvarrank*powerrank;
                realcorrs = reshape(tmp,wp.numfrex,nTimepoints)';
                PV(ian).CorrVals(nti,:,:,iv)  = realcorrs;
                PV(ian).time = linspace(-Pp.pwin(1),Pp.pwin(2),nTimepoints);
                
                
                % initialize null hypothesis matrices
                permuted_rvals  = zeros(wp.n_permutes,wp.numfrex,nTimepoints);
                max_pixel_rvals = zeros(wp.n_permutes,2);
                max_clust_info  = zeros(wp.n_permutes,1);

                % generate pixel-specific null hypothesis parameter distributions
                parfor permi = 1:wp.n_permutes
                    fake_rt_mapping = expvarrank(randperm(numrips));
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
                fprintf('running %d permutes\n', wp.n_permutes);
                tic
                parfor permi = 1:wp.n_permutes % i should par for this
                    % indices of permutations to include in thresholding at this iteration
                    perms2use4distribution = true(1,wp.n_permutes);
                    perms2use4distribution(permi) = 0;
                    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
                    fakecorrsz = squeeze((permuted_rvals(permi,:,:)-mean(permuted_rvals(perms2use4distribution,:,:),1)) ./ std(permuted_rvals(perms2use4distribution,:,:),[],1) );
                    fakecorrsz(abs(fakecorrsz)<norminv(1-wp.voxel_pval))=0;
                    % get number of elements in largest supra-threshold cluster
                    clustinfo = bwconncomp(fakecorrsz);
                    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
                end
                fprintf('took %d sec\n', toc);
                
                % now compute Z-map
                a = squeeze(PV(ian).CorrVals(nti,:,:,iv))';
                b = squeeze(mean(permuted_rvals,1));
                zmap = (a-b)./squeeze(std(permuted_rvals));
                PV(ian).zmap(nti,:,:,iv) = zmap;

                % apply pixel-level corrected threshold
                lower_threshold = prctile(max_pixel_rvals(:,1),wp.mcc_voxel_pval*100/2);
                upper_threshold = prctile(max_pixel_rvals(:,2),100-wp.mcc_voxel_pval*100/2);
                zmapthresh = zmap;
                zmapthresh(realcorrs>lower_threshold & realcorrs<upper_threshold)=0;
                PV(ian).pixelZmapThresh(nti,:,:,iv) = zmapthresh;

                % apply cluster-level corrected threshold
                zmapthresh = zmap;
                % uncorrected pixel-level threshold
                zmapthresh(abs(zmapthresh)<norminv(1-wp.voxel_pval))=0;
                % find islands and remove those smaller than cluster size threshold
                clustinfo = bwconncomp(zmapthresh);
                clust_info = cellfun(@numel,clustinfo.PixelIdxList);
                clust_threshold = prctile(max_clust_info,100-wp.mcc_cluster_pval*100);
                % identify clusters to remove
                whichclusters2remove = find(clust_info<clust_threshold);
                % remove clusters
                for i=1:length(whichclusters2remove)
                    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
                end
                PV(ian).clusterZmapThresh(nti,:,:,iv) = zmapthresh;
        end
    end
    end
end
if save_regression
    save_data(PV, [pconf.andef{2}, 'powerVarCorr/'], 'powerVarCorr_wtrack');
end
fprintf('started 9:55pm, ended %s\n', datetime)
if load_regression
    load_data([pconf.andef{2}, 'powerVarCorr/'], 'powerVarCorr_wtrack', Fp.animals);
end
if plot_regression
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({pwr.animal}, animal));
        for iv = 1:length(designmat(ian).expvars)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            set(gcf,'color','white')
            numrows = 2;
            numcols = ceil(length(ntrodes) / 2);
            
            for nti = 1:length(ntrodes)
                sf = subaxis(numrows,numcols,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    Pp.MgTp, 'MarginBottom', Pp.MgBm);
                nt = ntrodes(nti);
                area = ntinfo{1}{1}{nt}.area;
                subarea = ntinfo{1}{1}{nt}.subarea;
                if isnumeric(subarea)
                    subarea = num2str(subarea);
                end
                
                zmap = PV(ian).zmap(nti, :,:,iv);
                zmapthresh = PV(ian).clusterZmapThresh(nti, :,:,iv);
                
                    
                contourf(PV(ian).time,PV(ian).frequency,zmap,40,'linecolor','none')
                hold on
                 [~,h] = contour(PV(ian).time,PV(ian).frequency,logical(zmapthresh),1)            
                  h.LineColor = 'black';
                set(gca,'ydir','normal','yscale','log');
                caxis(sf, 'auto')
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end            
                            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s %s ', animal, Fp.uselfptype,designmat(ian).expvars{iv});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                pconf = animaldef('Demetris');
                save_figure(sprintf('%s/powerVarCorr/',pconf{4}), 'powerVarCorr', sprtit)
                close all
            end
            close all;
            %                 end
        end
    end
                
end

%% get the r2 and pval for TFboxVar for each nt, var
if run_TFVarcorr
    % load the raw power and do a tiedrank of the largebox meanpower across trials
        Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));

        % TF box
        PV(ian).TFboxidx = [fmnIdx, fmxIdx, tmnIdx, tmxIdx]; %bottom top left right
        PV(ian).Fminmax = [wp.fmn wp.fmx]; % user input vars
        PV(ian).Tminmax = [wp.tmn wp.tmx]; % user input vars
        % convert from Hz, sec to knn indices based on time, frex vecs
%         fmnIdx = 
%         fmxIdx =
%         tmnIdx =
%         tmxIdx = 

        for iv = 1:length(designmat(ian).expvars)
            fprintf('expvar %d\n', iv);
            expvarrank = designmat(ia).dm(:,iv)'; 
            numrips = length(find(~isnan(expvarrank)));
            expvarrank = expvarrank(~isnan(expvarrank));
            include_rips = ripstate.statesets(:,1);
            include_rips = find(~isnan(expvarrank));
            for nti = 1:length(ntrodes)
                fprintf('nt %d ',nti);
                
                % get mean pwr within tfbox for each rip for this nt
                TFBoxMeanRip = mean(mean(rawpwr(ian).data(in,tmnIdx:tmxIdx,:,fmnIdx:fmxIdx)));
                PV(ian).TFBoxCC(in,iv) = corrcoeff(TFBoxMeanRip, designmat(ian).data(:,iv));

                %  screw this small box thing.. once i've picked the coords for the large
                %  box.. run corrcoef for the mean rawpwr within the large box vs var

        %         % find peak within TFLargeidx for this nt iv tf 
        %         largeBoxData = squeeze(PV(ian).data(in,tmnIdx:tmxIdx,fmnIdx:fmxIdx,iv));
        %         [Mtf, I] = max(largeBoxData(:));
        %         [tcentIdx, fcentIdx] = ind2sub(Mtf, I);
        %         PV(ian).TpeakIdx = tmnIdx + tcentIdx;
        %         PV(ian).FpeakIdx = fmnIdx + fcentIdx;

        %         % small box size params 
        %         PV(ian).TFSmallidx = wp.TFSmallidx; % [10 10 2 2] bottom top left right
        %         
        %         smFmnidx = PV(ian).FpeakIdx - wp.TFSmallidx(1);
        %         smFmxidx = PV(ian).FpeakIdx + wp.TFSmallidx(2);
        %         smTmnidx = PV(ian).TpeakIdx - wp.TFSmallidx(3);
        %         smTmxidx = PV(ian).TpeakIdx + wp.TFSmallidx(4);

        %         PV(ian).smallTFcorr = mean(mean(PV(ian).data(in,smTmnidx:smFmxidx,smFmnidx:smFmxidx,iv));

        %         % get indices of small box centered at peak
        %         PV(ian).Fwidth = wp.frexSmallBoxWidthHz
        %         
        %         % convert from Hz, sec to knn indi
        %         
        %         % run corrcoef on tfpwr vec and var vec for each nt and each var
        %         PV(ian).pval(in,iv) = 
        %         PV(ian).r2(in,iv) = 
    end
    end
    end
end

if plot_corrheatmap
   % plot each animal's PV(ian).TFBoxCC heatmap for ntrode x vars
   % identify the nt's sig across anims for each var
   % is there a var that has sig nt for each an for the chosen tf?
   
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