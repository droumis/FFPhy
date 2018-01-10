function [baseLinemeanpower, pctbasedpowerout,logbasedpowerout] = computePower(as,ian,int,iEpTetIndsType,baseind, indsSet)
% [powerout, basedmeanpower, percbasedpowerout,logbasedpowerout, logbasedpoweroutUnspec,percbasedpoweroutUnspec,basedmeanpowerUnspec]  =...
%     computePower(as,ian,int,iEpTetIndsType,baseind, indsSet)
disp(sprintf('$$$$$$$$$$$$$$$ POWER for animal %s nt %s  $$$$$$$$$$$$$$$$$$$',ian, int));
% 
% % condition-unspecific baseline mean
% %     basedmeanpower = mean(mean(abs(ias).^2,2),1); %mean for each freq across time and across all events for the baseline period
%     basedmeanpowerUnspec = mean(mean(abs(as(baseind(1):baseind(2),:,:)).^2,2),1);
    
for iSet = 1:size(iEpTetIndsType, 2)
    
    
    %% compute power as the magnitude from origin in complex space
    % un-normalized power data
    powertmp = cellfun(@(x) mean(abs(as(:,x,:)).^2,2), iEpTetIndsType(:,iSet), 'un', 0); %get the mean across events for each event state type
    powerout{iSet} = cat(2,powertmp{:}); %col cat the means for each event state type
    
    %% baseline data, mean
    % baselineASignal = as(baseind(1):baseind(2),iEpTetIndsType{1},:); %analytic signal from baseline period
   
    
    % condition-specific baseline mean
    powerbaselinemeans = cellfun(@(x) mean(mean(abs(as(baseind(1):baseind(2),x,:)).^2,2),1),iEpTetIndsType(:,iSet), 'un', 0);
    % powerbaselinemeans = cellfun(@(x) mean(mean(abs(baselineASignal(:,x,:)).^2,2),1),iEpTetIndsType, 'un', 0);
    baseLinemeanpower{iSet} = cat(2,powerbaselinemeans{:});
    
    %% Normalization
    
    % % z-score normalization using only the baseline window to compute mean, std; condition-specific baseline mean
    % powerbaselinestds = cellfun(@(x) std(reshape(abs(ias(:,x,:)).^2,[],1,numfrex),1),iEpTetIndsType, 'un', 0);
    % basedstdpower = cat(2,powerbaselinestds{:});
    % zscoredbasedpowerOut = bsxfun(@rdivide, bsxfun(@minus, powerOut,basedmeanpower), basedstdpower);
    %
    % % z-score normalization using the entire window to compute mean, std;  condition-specific baseline mean
    % % vert cat all the traces into columns of frex, then std within freq of the baseline period
    % powermeans = cellfun(@(x) mean(mean(abs(as{int}(:,x,:)).^2,2),1),iEpTetIndsType, 'un', 0);
    % meanpower = cat(2,powermeans{:});
    % powerstds = cellfun(@(x) std(reshape(abs(as{int}(:,x,:)).^2,[],1,numfrex),1),iEpTetIndsType, 'un', 0);
    % stdpower = cat(2,powerstds{:});
    % zscoredpowerOut = bsxfun(@rdivide, bsxfun(@minus, powerOut,meanpower), stdpower);
    %
    % % percent change from baseline - normalization using condition-specific baseline mean
    pctbasedpowerout{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, powerout{iSet},baseLinemeanpower{iSet}), baseLinemeanpower{iSet});
    
    % decibel - normalization using condition-specific baseline mean
    logbasedpowerout{iSet} = 10*log10(bsxfun(@rdivide, powerout{iSet},baseLinemeanpower{iSet}));
    
%     % % percent change from baseline - normalization using condition-UNspecific baseline mean
%     percbasedpoweroutUnspec{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, powerout{iSet},basedmeanpowerUnspec), basedmeanpowerUnspec);
%     
%     % decibel - normalization using condition-UNspecific baseline mean
%     logbasedpoweroutUnspec{iSet} = 10*log10(bsxfun(@rdivide, powerout{iSet},basedmeanpowerUnspec));
%     
    %% compute differential tf maps
    if strcmp(indsSet, 'performance') || strcmp(indsSet, 'performanceByDay')
        try
            powerout{iSet}(:,6,:) = powerout{iSet}(:,4,:) - powerout{iSet}(:,5,:);
            powerout{iSet}(:,7,:) = powerout{iSet}(:,2,:) - powerout{iSet}(:,3,:);
            logbasedpowerout{iSet}(:,6,:) = logbasedpowerout{iSet}(:,4,:) - logbasedpowerout{iSet}(:,5,:);
            logbasedpowerout{iSet}(:,7,:) = logbasedpowerout{iSet}(:,2,:) - logbasedpowerout{iSet}(:,3,:);
            pctbasedpowerout{iSet}(:,6,:) = pctbasedpowerout{iSet}(:,4,:) - pctbasedpowerout{iSet}(:,5,:);
            pctbasedpowerout{iSet}(:,7,:) = pctbasedpowerout{iSet}(:,2,:) - pctbasedpowerout{iSet}(:,3,:);
            % % i actually don't think you want to use the raw difference. instead use
            % the already normalized, then differenced
            % tmpinds = [iEpTetIndsType{4}; iEpTetIndsType{5}];
            % basedmeanpower(:,6,:) = mean(mean(abs(baselineASignal(:,tmpinds,:)).^2,2),1);
            % tmpinds = [iEpTetIndsType{2}; iEpTetIndsType{3}];
            % basedmeanpower(:,7,:) = mean(mean(abs(baselineASignal(:,tmpinds,:)).^2,2),1);
            
%             logbasedpoweroutUnspec{iSet}(:,6,:) = logbasedpoweroutUnspec{iSet}(:,4,:) - logbasedpoweroutUnspec{iSet}(:,5,:);
%             logbasedpoweroutUnspec{iSet}(:,7,:) = logbasedpoweroutUnspec{iSet}(:,2,:) - logbasedpoweroutUnspec{iSet}(:,3,:);
%             percbasedpoweroutUnspec{iSet}(:,6,:) = percbasedpoweroutUnspec{iSet}(:,4,:) - percbasedpoweroutUnspec{iSet}(:,5,:);
%             percbasedpoweroutUnspec{iSet}(:,7,:) = percbasedpoweroutUnspec{iSet}(:,2,:) - percbasedpoweroutUnspec{iSet}(:,3,:);
% 
        catch
            disp(' ..')
        end
    end
end

end

















































































