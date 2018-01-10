function[ISPCout, basedmeanISPC, pctbasedISPCout, logbasedISPCout] = computeISPC(ph,ian,int,iEpTetIndsType,baseind, indsSet, varargin)

wp = [];
if ~isempty(varargin)
    assign(varargin{:});
end
% [ISPCout, basedmeanISPC, pctbasedISPCout, logbasedISPCout, logbasedISPCoutUnspec, pctbasedISPCoutUnspec, unspecBaselineISPCmean] =...
%     computeISPC(ph,ian,int,iEpTetIndsType,baseind, indsSet)

disp(sprintf('$$$$$$$$$$$$$$$ ISPC for animal %s nt%d  $$$$$$$$$$$$$$$$$$$',ian, int));

% condition-unspecific baseline mean
% unspecBaselineISPCmean = mean(abs(mean(exp(1i*ph(baseind(1):baseind(2),:,:)),2)),1);

% initialize output time-frequency data
ixpc.output{ianimal}{iday} = zeros(nchoosek(size(iEpTetDataCat,1),2), wp.nsamps, wp.num_frex);
ixpc.phOut{ianimal}{idatatype} = zeros(nchoosek(size(iEpTetDataCat,1),2), wp.nsamps, wp.num_frex);



for iSet = 1:size(iEpTetIndsType, 1)
    % if strcmp(indsSet, 'performanceByDay')
    
    
    %% compute power as the magnitude from origin in complex space
    % un-normalized phase data
    ISPCtmp = cellfun(@(x) abs(mean(exp(1i*ph(:,x,:)),2)), iEpTetIndsType(iSet, :), 'un', 0); %get the mean across events for each event state type
    %     ixpc.ISPCout{ian}{int}(:,:,:) = cat(2, ISPCtmp{:}); %col cat the means for each event state type
    ISPCout{iSet} = cat(2, ISPCtmp{:}); %col cat the means for each event state type
    %% baseline data, mean
    % baselinePhData = cellfun(@(x) exp(1i*ph(baseind(1):baseind(2),iEpTetIndsType{x},:)), baselinePhDataIndSet, 'un', 0);
    
    % condition-specific baseline mean
    baselinesISPCmean = cellfun(@(x) mean(abs(mean(exp(1i*ph(baseind(1):baseind(2),x,:)),2)),1),iEpTetIndsType(iSet, :), 'un', 0);
    
    % baselinesISPCmean = cellfun(@(x) mean(abs(mean(baselinePhData(:,x,:),2)),1),iEpTetIndsType, 'un', 0);
    %     ixpc.basedmeanISPC{ian}{int} = cat(2,baselinesISPCmean{:});
    basedmeanISPC{iSet} = cat(2,baselinesISPCmean{:});
    
    %% Normalization
    
    % % z-score normalization using only the baseline window to compute mean, std; condition-specific baseline mean
    
    % % z-score normalization using the entire window to compute mean, std;  condition-specific baseline mean
    % % vert cat all the traces into columns of frex, then std within freq of the baseline period
    
    % % percent change from baseline - normalization using condition-specific baseline mean
    %     ixpc.percbasedISPCout{ian}{int} = 100 * bsxfun(@rdivide, bsxfun(@minus, ixpc.ISPCout{ian}{int},ixpc.basedmeanISPC{ian}{int}), ixpc.basedmeanISPC{ian}{int});
    pctbasedISPCout{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, ISPCout{iSet},basedmeanISPC{iSet}), basedmeanISPC{iSet});
    
    % percent change from baseline - normalization using condition-UNspecific baseline mean
%     pctbasedISPCoutUnspec{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, ISPCout{iSet},unspecBaselineISPCmean), basedmeanISPC{iSet});
    
    % decibel - normalization using condition-specific baseline mean
    %     ixpc.logbasedISPCout{ian}{int} = 10*log10(bsxfun(@rdivide, ixpc.ISPCout{ian}{int},ixpc.basedmeanISPC{ian}{int}));
    logbasedISPCout{iSet} = 10*log10(bsxfun(@rdivide, ISPCout{iSet},basedmeanISPC{iSet}));
    
    % decibel - normalization using condition-UNspecific baseline mean
%     logbasedISPCoutUnspec{iSet} = 10*log10(bsxfun(@rdivide, ISPCout{iSet},unspecBaselineISPCmean));
    
    % ixpc.basedstdphase{ian}{int}(1,1,:) = std(reshape(abs(mean(iph,2)),[],1,numfrex),1);
    
    % baseline Z- normalize. don't use this anymore
    % ixpc.basedphOut{ian}{int} = bsxfun(@rdivide, bsxfun(@minus, ixpc.phOut{ian}{int},ixpc.basedmeanphase{ian}{int}), ixpc.basedstdphase{ian}{int});
    % percent change from baseline
    
    % ixpc.zscoredphOut{ian}{int} = bsxfun(@rdivide, bsxfun(@minus, ixpc.phOut{ian}{int},ixpc.basedmeanphase{ian}{int}), ixpc.basedmeanphase{ian}{int});
    
    %% compute differential tf maps
    if strcmp(indsSet, 'performance') || strcmp(indsSet, 'performanceByDay')
        try
            ISPCout{iSet}(:,6,:) = ISPCout{iSet}(:,4,:) - ISPCout{iSet}(:,5,:);
            ISPCout{iSet}(:,7,:) = ISPCout{iSet}(:,2,:) - ISPCout{iSet}(:,3,:);
            % i actually don't think you want to use the raw difference. instead use
            % the already normalized, then differenced
            logbasedISPCout{iSet}(:,6,:) = logbasedISPCout{iSet}(:,4,:) - logbasedISPCout{iSet}(:,5,:);
            logbasedISPCout{iSet}(:,7,:) = logbasedISPCout{iSet}(:,2,:) - logbasedISPCout{iSet}(:,3,:);
            pctbasedISPCout{iSet}(:,6,:) = pctbasedISPCout{iSet}(:,4,:) - pctbasedISPCout{iSet}(:,5,:);
            pctbasedISPCout{iSet}(:,7,:) = pctbasedISPCout{iSet}(:,2,:) - pctbasedISPCout{iSet}(:,3,:);
    
%             % condition unspecific baseline normalized diff
%             logbasedISPCoutUnspec{iSet}(:,6,:) = logbasedISPCoutUnspec{iSet}(:,4,:) - logbasedISPCoutUnspec{iSet}(:,5,:);
%             logbasedISPCoutUnspec{iSet}(:,7,:) = logbasedISPCoutUnspec{iSet}(:,2,:) - logbasedISPCoutUnspec{iSet}(:,3,:);
%             pctbasedISPCoutUnspec{iSet}(:,6,:) = pctbasedISPCoutUnspec{iSet}(:,4,:) - pctbasedISPCoutUnspec{iSet}(:,5,:);
%             pctbasedISPCoutUnspec{iSet}(:,7,:) = pctbasedISPCoutUnspec{iSet}(:,2,:) - pctbasedISPCoutUnspec{iSet}(:,3,:);
        catch
            disp(' ..')
        end
    end
end
end




    
    
    %                 phdata{introde}(:,:,fi) = angle(as{introde}(:,:,fi)); %get phase component of the analytic signal
    %             %% ISPC -- NEED TO FIX :(
    %             if calculateISPC % time series diff of all possible ntrode pairs.
    %                 %the idea here is to treat rows as pairs instead of individual ntrodes, as in ITPC
    %                 ispc = zeros(nchoosek(size(phdata{introde}{fi},1),2), nsamps, nevents);
    %                                 m = logical(tril(ones(size(phdata{introde}{fi}(:,:,1),1)),-1)); %get indices of non-duplicates (below comb triangle)
    %                                 [brow, bcol] = find(m); %get linear index of ntrode combination indices
    %                                 indices = [F(ianimal).output{day}(eps(1)).index(bcol,:) F(ianimal).output{day}(eps(1)).index(brow,:)]; %convert to ntrodeID
    %                 for w = 1:size(phdata{introde}{fi},3); %loop through each event plane (ntrode x samples X event#)
    %                     iprm = permute(phdata{introde}{fi}(:,:,w),[3 2 1]); %transpose
    %                     B = bsxfun(@minus,phdata{introde}{fi}(:,:,w),iprm); %subtract across each NTrode-choose-two comb of phase time series
    %                     B = reshape(B,[],size(phdata{introde}{fi}(:,:,w),2)); %reshape back into 2D
    %                     B = B(m(:),:); %get rid of duplicates
    %                     ispc(:,:,w) = B; %save result
    %                 end
    %                 clear phasedata
    %                 phdata{introde}{fi} = ispc;
    %             end
    %             ixpc.index{ianimal}{fi} = indices;
    


















































































