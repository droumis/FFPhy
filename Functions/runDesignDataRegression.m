

function PV = runDesignDataRegression(design, rawpwr, Fp, varargin)
% inputs:
% design struct containing design matrix
% rawpwr struct with dims ntrode, time, rips, frex

% outputs:
% realbeta for each nt, time, frex, expvar
% zmap ""
% zmapthresh ""
invalid_rips = struct;
pconf = paramconfig;
saveout = 1;
outdir = 'varCont';
if ~isempty(varargin)
    assign(varargin{:})
end
outpath = [pconf.andef{2},outdir,'/'];
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
    pwranidx = find(strcmp({rawpwr.animal}, animal));
    if ~isempty(invalid_rips)
        % exclude invalid rips
        noiseanidx = find(strcmp({invalid_rips.animal}, animal));
        invalidrips = invalid_rips(noiseanidx ).ripnums;
        userips = ones(length(rawpwr(pwranidx).day),1);
        userips(invalidrips) = 0;
    end
    
    PV(ian).animal = animal;
    PV(ian).designMat = design(desani);
    PV(ian).ntrodes = ntrodes;
    PV(ian).expvars = design(desani).expvars;
    PV(ian).frequency = wp.frex;
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
        % get in,iv tfcorr for design vs data
        %         [PV(ian).data(in,:,:,iv), PV(ian).zmap(in,:,:,iv), PV(ian).thresh(in,:,:,:)] = ...
        %             regress(designmat(ian).data(:,iv), rawpwr(ian).data(nt,win1:win2,:));
        % rotate to make time dim 2
        data = permute(squeeze(rawpwr(pwranidx).pwr(nti, :, include_rips,:)),[2 1 3]);
        data = trim2win(data, Fp.srate, Pp.pwin, 'dsamp', wp.dsamp);
        nTimepoints = size(data,2);
        data = permute(data, [3 2 1]); % rotate back to [freq time ripple]
        powerreshaped = reshape(data, wp.numfrex*nTimepoints,numrips); % flatten [ripple timefreq]
        powerrank = tiedrank(powerreshaped'); % tiedrank on dim 1, across rips.
        % powerrank (1168, 37525) numbers 1:1168
        % X (25, 1168) numbers -2:16, mostly -2:4
%         a = (X*X'); %not nan (25, 25) numbers tens to hundreds. thousand on diag
%         b = powerrank; % b is not nan (25, 37525) numbers in the thousands, tens thousands
%         realbeta = a\b; % EVERYTHING IS NAN.. are the b numbers too large as denominators?
        % something is going on with the line above to make everything Nan... 
%         realbeta = (X*X')\X; % aha!! this is nan.. so it's a design matrix issue
        % aha! the issue was that i was using too many tfboxes.. since this
        % is multiregression, it was making the beta's infinitely small
        realbeta = (X*X')\X*powerrank;
        realbeta = reshape(realbeta,numvars,wp.numfrex,nTimepoints);
        PV(ian).CorrVals(nti,:,:,:)  = realbeta;
        PV(ian).time = linspace(-Pp.pwin(1),Pp.pwin(2),nTimepoints);
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
%             a = fakeX*powerrank;
%             b = fakebeta\a;
%             if sum(sum(isnan(b)))
%                 fprintf('nans in the pudding\n')
% %                 pause
%             end
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
            PV(ian).zmap(nti,:,:,testi) = zmap;
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
            PV(ian).clusterZmapThresh(nti,:,:,testi) = zmapthresh;
        end
    end
    
    if saveout
        save_data(PV(ian), outpath, [outdir,'Corr_',Fp.epochEnvironment]);
    end
end
fprintf('took %d sec\n', toc);
end