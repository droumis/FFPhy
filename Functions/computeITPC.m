function[baseLineITPC, pctbasedITPCout, logbasedITPCout] = computeIXPC(ph,ian,int,iEpTetIndsType,baseind, indsSet)
% [ITPCout, basedmeanITPC, pctbasedITPCout, logbasedITPCout, logbasedITPCoutUnspec, pctbasedITPCoutUnspec, unspecBaselineITPCmean] =...
%     computeITPC(ph,ian,int,iEpTetIndsType,baseind, indsSet)

disp(sprintf('$$$$$$$$$$$$$$$ ITPC for animal %s nt%d  $$$$$$$$$$$$$$$$$$$',ian, int));

% condition-unspecific baseline mean
% unspecBaselineITPCmean = mean(abs(mean(exp(1i*ph(baseind(1):baseind(2),:,:)),2)),1);
    
for iSet = 1:size(iEpTetIndsType, 2)
    % if strcmp(indsSet, 'performanceByDay')
    
    
    %% compute power as the magnitude from origin in complex space
    % un-normalized phase data
    itpctmp = cellfun(@(x) abs(mean(exp(1i*ph(:,x,:)),2)), iEpTetIndsType(:, iSet), 'un', 0); %get the mean across events for each event state type
    %     ixpc.ITPCout{ian}{int}(:,:,:) = cat(2, itpctmp{:}); %col cat the means for each event state type
    ITPCout{iSet} = cat(2, itpctmp{:}); %col cat the means for each event state type
    %% baseline data, mean
    % baselinePhData = cellfun(@(x) exp(1i*ph(baseind(1):baseind(2),iEpTetIndsType{x},:)), baselinePhDataIndSet, 'un', 0);
    
    % condition-specific baseline mean
    baselinesITPCmean = cellfun(@(x) mean(abs(mean(exp(1i*ph(baseind(1):baseind(2),x,:)),2)),1),iEpTetIndsType(:, iSet), 'un', 0);
    
    % baselinesITPCmean = cellfun(@(x) mean(abs(mean(baselinePhData(:,x,:),2)),1),iEpTetIndsType, 'un', 0);
    %     ixpc.basedmeanITPC{ian}{int} = cat(2,baselinesITPCmean{:});
    baseLineITPC{iSet} = cat(2,baselinesITPCmean{:});
    
    %% Normalization
    
    % % z-score normalization using only the baseline window to compute mean, std; condition-specific baseline mean
    
    % % z-score normalization using the entire window to compute mean, std;  condition-specific baseline mean
    % % vert cat all the traces into columns of frex, then std within freq of the baseline period
    
    % % percent change from baseline - normalization using condition-specific baseline mean
    %     ixpc.percbasedITPCout{ian}{int} = 100 * bsxfun(@rdivide, bsxfun(@minus, ixpc.ITPCout{ian}{int},ixpc.basedmeanITPC{ian}{int}), ixpc.basedmeanITPC{ian}{int});
    pctbasedITPCout{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, ITPCout{iSet},baseLineITPC{iSet}), baseLineITPC{iSet});
    
    % percent change from baseline - normalization using condition-UNspecific baseline mean
%     pctbasedITPCoutUnspec{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, ITPCout{iSet},unspecBaselineITPCmean), basedmeanITPC{iSet});
    
    % decibel - normalization using condition-specific baseline mean
    %     ixpc.logbasedITPCout{ian}{int} = 10*log10(bsxfun(@rdivide, ixpc.ITPCout{ian}{int},ixpc.basedmeanITPC{ian}{int}));
    logbasedITPCout{iSet} = 10*log10(bsxfun(@rdivide, ITPCout{iSet},baseLineITPC{iSet}));
    
    % decibel - normalization using condition-UNspecific baseline mean
%     logbasedITPCoutUnspec{iSet} = 10*log10(bsxfun(@rdivide, ITPCout{iSet},unspecBaselineITPCmean));
    
    % ixpc.basedstdphase{ian}{int}(1,1,:) = std(reshape(abs(mean(iph,2)),[],1,numfrex),1);
    
    % baseline Z- normalize. don't use this anymore
    % ixpc.basedphOut{ian}{int} = bsxfun(@rdivide, bsxfun(@minus, ixpc.phOut{ian}{int},ixpc.basedmeanphase{ian}{int}), ixpc.basedstdphase{ian}{int});
    % percent change from baseline
    
    % ixpc.zscoredphOut{ian}{int} = bsxfun(@rdivide, bsxfun(@minus, ixpc.phOut{ian}{int},ixpc.basedmeanphase{ian}{int}), ixpc.basedmeanphase{ian}{int});
    
    %% compute differential tf maps
    if strcmp(indsSet, 'performance') || strcmp(indsSet, 'performanceByDay')
        try
            ITPCout{iSet}(:,6,:) = ITPCout{iSet}(:,4,:) - ITPCout{iSet}(:,5,:);
            ITPCout{iSet}(:,7,:) = ITPCout{iSet}(:,2,:) - ITPCout{iSet}(:,3,:);
            % i actually don't think you want to use the raw difference. instead use
            % the already normalized, then differenced
            logbasedITPCout{iSet}(:,6,:) = logbasedITPCout{iSet}(:,4,:) - logbasedITPCout{iSet}(:,5,:);
            logbasedITPCout{iSet}(:,7,:) = logbasedITPCout{iSet}(:,2,:) - logbasedITPCout{iSet}(:,3,:);
            pctbasedITPCout{iSet}(:,6,:) = pctbasedITPCout{iSet}(:,4,:) - pctbasedITPCout{iSet}(:,5,:);
            pctbasedITPCout{iSet}(:,7,:) = pctbasedITPCout{iSet}(:,2,:) - pctbasedITPCout{iSet}(:,3,:);
    
%             % condition unspecific baseline normalized diff
%             logbasedITPCoutUnspec{iSet}(:,6,:) = logbasedITPCoutUnspec{iSet}(:,4,:) - logbasedITPCoutUnspec{iSet}(:,5,:);
%             logbasedITPCoutUnspec{iSet}(:,7,:) = logbasedITPCoutUnspec{iSet}(:,2,:) - logbasedITPCoutUnspec{iSet}(:,3,:);
%             pctbasedITPCoutUnspec{iSet}(:,6,:) = pctbasedITPCoutUnspec{iSet}(:,4,:) - pctbasedITPCoutUnspec{iSet}(:,5,:);
%             pctbasedITPCoutUnspec{iSet}(:,7,:) = pctbasedITPCoutUnspec{iSet}(:,2,:) - pctbasedITPCoutUnspec{iSet}(:,3,:);
        catch
            disp(' ..')
        end
    end
end
end






















































































