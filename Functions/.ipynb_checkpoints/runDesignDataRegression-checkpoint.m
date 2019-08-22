

function out = runDesignDataRegression(design, rawpwr, Fp, varargin)
% inputs:
% design struct containing design matrix
% rawpwr struct with dims ntrode, time, rips, frex
% outputs:
% realbeta for each nt, time, frex, expvar
% zmap ""
% zmapthresh ""
noiseEvents = struct;
pconf = paramconfig;
saveout = 1;
outdir = 'varCont';
tfboxes = []; 
if ~isempty(varargin)
    assign(varargin{:})
end
tic
Pp = load_plotting_params({'defaults', 'power'});
wp = getWaveParams('4-300Hz');
fprintf('running multi regression with %d permutes\n', wp.n_permutes);
for ian = 1:length(rawpwr)
    animal = rawpwr(ian).animal;
    desani = find(strcmp({design.animal}, animal));
    aninfo = animaldef(animal);
    ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
    ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
    ntrodes = unique(ntrodes(:,3));
    userips = ones(length(rawpwr(ian).day),1);
    if ~isempty(noiseEvents)
        % exclude invalid rips
        noiseanidx = find(strcmp({noiseEvents.animal}, animal));
        if ~isempty(noiseEvents(noiseanidx).events)
            invalidrips = ismember([rawpwr(ian).day rawpwr(ian).epoch rawpwr(ian).ripStartTime], ...
                noiseEvents(noiseanidx).events, 'rows');
            userips(invalidrips) = 0;
        end
    end
    out(ian).animal = animal;
    out(ian).dm = design(desani);
    out(ian).ntrodes = ntrodes;
    out(ian).expvars = design(desani).expvars;
    out(ian).frequency = wp.frex;
    %         for iv = 1:length(designmat(ian).expvars)
    %             fprintf('expvar %d\n', iv);
    % only keep the rips that have non nan vals in all the vars
    numvars = length(design(desani).expvars);
    if length(size(design(desani).dm)) == 2
        include_rips = all(~isnan(design(desani).dm),2);
        include_rips(~userips) = 0;
        numrips = sum(include_rips(:,1));
        X = zscore(design(desani).dm(include_rips,:))';
        numrips = sum(include_rips(:,1));
        numvars = length(X(:,1));
    end
    %         expvarZ = expvarZ(~isnan(expvarZ));
    %             include_rips = ones(length(expvarrank),1);
    %         include_rips = find(~isnan(expvarZ));
    
    for nti = 1:length(ntrodes)
        fprintf('nt %d \n',nti);
        if length(size(design(desani).dm)) == 3
            include_rips = squeeze(all(all(~isnan(design(desani).dm),2),3));
            include_rips(~userips) = 0;
            numrips = sum(include_rips(:,1));
            X = zscore(squeeze(design(desani).dm(include_rips,:,nti)))';
        end
        % rotate to make time dim 2
        data = permute(squeeze(rawpwr(ian).pwr(nti, :, include_rips,:)),[2 1 3]);
        data = trim2win(data, Fp.srate, Pp.pwin, 'dsamp', wp.dsamp);
        nTimepoints = size(data,2);
        data = permute(data, [3 2 1]); % rotate back to [freq time ripple]
        powerreshaped = reshape(data, wp.numfrex*nTimepoints,numrips); % flatten [ripple timefreq]
        powerrank = tiedrank(powerreshaped'); % tiedrank on dim 1, across rips.
        realbeta = (X*X')\X*powerrank; % regression
        realbeta = reshape(realbeta,numvars,wp.numfrex,nTimepoints);
        out(ian).CorrVals(nti,:,:,:)  = realbeta;
        out(ian).time = linspace(-Pp.pwin(1),Pp.pwin(2),nTimepoints);
        % initialize null hypothesis matrices
        permuted_bvals  = zeros(wp.n_permutes,numvars,wp.numfrex,nTimepoints);
        %             max_pixel_rvals = zeros(wp.n_permutes,numvars);
        max_clust_info  = zeros(wp.n_permutes,numvars);
        % generate pixel-specific null hypothesis parameter distributions
        for permi = 1:wp.n_permutes %parfor this
            % randomly shuffle trial order
            fakeX = X(:,randperm(numrips));
            % compute beta-map of null hypothesis.. i.e. do the
            % multiregression on the ripple-shuffled design matrix and power
            fakebeta = (fakeX*fakeX')\fakeX*powerrank;
            % reshape to 2D map for cluster-correction
            fakebeta = reshape(fakebeta,numvars,wp.numfrex,nTimepoints);
            % save all permuted values
            permuted_bvals(permi,:,:,:) = fakebeta;
            % save maximum pixel values
            %                 max_pixel_rvals(permi,:) = [ min(fakebeta(:)) max(fakebeta(:)) ];
        end
        % the cluster correction will be done on the permuted data, thus
        % making no assumptions about parameters for p-values
        for permi = 1:wp.n_permutes % i should par for this
            for testi = 1:numvars
                % for cluster correction, apply uncorrected threshold and
                % get maximum cluster sizes
                fakecorrsz = squeeze((permuted_bvals(permi,testi,:,:)- ...
                    mean(permuted_bvals(:,testi,:,:),1)) ./ ...
                    std(permuted_bvals(:,testi,:,:),[],1) );
                fakecorrsz(abs(fakecorrsz)<norminv(1-wp.voxel_pval))=0;
                % get number of elements in largest supra-threshold cluster
                clustinfo = bwconncomp(fakecorrsz);
                % the zero accounts for empty maps
                max_clust_info(permi,testi) = max([ 0 ...
                    cellfun(@numel,clustinfo.PixelIdxList) ]);
            end
        end
        for testi = 1:numvars % parfor this
            % now compute Z-map
            zmap = (squeeze(realbeta(testi,:,:))-...
                squeeze(mean(permuted_bvals(:,testi,:,:),1))) ./ ...
                squeeze(std(permuted_bvals(:,testi,:,:),[],1));
            out(ian).zmap(nti,:,:,testi) = zmap;
            % apply cluster-level corrected threshold
            zmapthresh = zmap;
            % uncorrected pixel-level threshold
            zmapthresh(abs(zmapthresh)<norminv(1-wp.voxel_pval))=0;
            % find islands and remove those smaller than cluster size threshold
            clustinfo = bwconncomp(zmapthresh);
            clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            clust_threshold = prctile(max_clust_info(:,testi),100-wp.mcc_cluster_pval*100);
            % identify clusters to remove
            whichclusters2remove = find(clust_info<clust_threshold);
            % remove clusters
            for i=1:length(whichclusters2remove)
                zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end
            out(ian).clusterZmapThresh(nti,:,:,testi) = zmapthresh;
        end
        tfbanim = strcmp({tfboxes.animal}, animal);
        if ~isempty(tfboxes)
            out(ian).tfboxes = tfboxes(tfbanim);
            %         % now compute correlation of the pwr mean of each tfbox and each
            %         % condition
            
            for ivar = 1:numvars
                if length(size(design(desani).dm)) == 3
                    Xvar = design(desani).dm(include_rips,ivar,nti)';
                else
                    Xvar = design(desani).dm(include_rips,ivar)';
                end
                for itfb = 1:length(tfboxes(tfbanim).expvars)
                    out(ian).tfb{itfb} = tfboxes(tfbanim).expvars{itfb};
                    itfbpwr = squeeze(tfboxes(tfbanim).dm(include_rips, itfb, nti))';
                    tmp = fitlm(Xvar, itfbpwr);
                    out(ian).fitlm{nti, itfb, ivar} = tmp;
                    out(ian).P(nti, itfb, ivar) = tmp.coefTest;
                    out(ian).R(nti, itfb, ivar) = tmp.Rsquared.Ordinary;
                    out(ian).coef(nti, itfb, ivar) = tmp.Coefficients.Estimate(end);
                end
            end
        end
    end
    
    if saveout
        outpath = [pconf.andef{2},outdir,'/'];
        save_data(out(ian), outpath, [outdir,'Corr_',Fp.epochEnvironment]);
    end
end
fprintf('took %d sec\n', toc);
end