
function ixpc = permutationTest(ixpc, as, ph, ian, int, IndicesByType, DataTypeFields, n_permutes, permTestType)

[condAInds, condBInds, outInd] = permTestParams(DataTypeFields, permTestType, int,n_permutes);

ixpc.powerzmask{ian}{int}{outInd} = []; ixpc.ITPCzmask{ian}{int}{outInd} = [];
ixpc.MC_power_minmax{ian}{int}{outInd} = []; ixpc.MC_ITPC_minmax{ian}{int}{outInd} = [];


[ixpc.powerzmask{ian}{int}{outInd}, ixpc.ITPCzmask{ian}{int}{outInd},...
    ixpc.MC_power_minmax{ian}{int}{outInd},ixpc.MC_ITPC_minmax{ian}{int}{outInd}]...
    = ITPCpermtest(ixpc.powerout{ian}{int}(:,outInd,:), ixpc.ITPCout{ian}{int}(:,outInd,:),...
    ph, as, IndicesByType{condAInds}, IndicesByType{condBInds}, n_permutes);
end


function [condAInds, condBInds, outInd] = permTestParams(DataTypeFields, permTestType, int, n_permutes)

switch permTestType
    %                 iEpTetDataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'};
    case 'outboundVSinbound'
        % outbound vs inbound
        disp(sprintf('========== running perm tests Outb vs Inb ============= nt%d perm x %d',int, n_permutes))
        condAInds = find(cell2mat(cellfun(@(x) strcmp(x,'outB'), DataTypeFields, 'UniformOutput', false)));
        condBInds = find(cell2mat(cellfun(@(x) strcmp(x,'inB'), DataTypeFields, 'UniformOutput', false)));
        outInd = find(cell2mat(cellfun(@(x) strcmp(x,'outB-inB'), DataTypeFields, 'UniformOutput', false)));
        
    case 'correctOutvsmistakeOut'
        % correct-outbound vs mistake-outbound
        disp(sprintf('========== running perm tests CorrectOutb vs MistakeOutb ============= nt%d perm x %d',int, n_permutes))
        condAInds = find(cell2mat(cellfun(@(x) strcmp(x,'corrOut'), DataTypeFields, 'UniformOutput', false)));
        condBInds = find(cell2mat(cellfun(@(x) strcmp(x,'mistOut'), DataTypeFields, 'UniformOutput', false)));
        outInd = find(cell2mat(cellfun(@(x) strcmp(x,'corrOut-mistOut'), DataTypeFields, 'UniformOutput', false)));
        
    otherwise
        error('pick an existing perm test to run')
        % ixpc = ITPCpermtest(ixpc, phdata{introde}, as{introde}, ianimal, introde, 7, iEpTetIndsType{2}, iEpTetIndsType{3}, n_permutes);
        %% cluster-based multiple comp scratch
        %                 for iprm = 1:n_permutes;
        % %                     % take each permutation map, and transform to Z
        % %                     ipermmap = phasepermOutput(:,iprm,:);
        % %                     izpermmap = (ipermmap-ixpc.phasemean_h0{introde}{6}(:,1,:))./ixpc.pwrstd_h0{introde}{6}(:,1,:);
        % %                     % threshold image at p-value
        % %                     izpermmap(abs(izpermmap)<zval) = 0;
        %                     % find clusters (need image processing toolbox for this!)
        %                     izpermmap = squeeze(zpermmaps(:,iprm,:));
        %                     izpermislands = bwconncomp(izpermmap);
        %                     if ~isempty(cell2mat(cellfun(@numel,izpermislands.PixelIdxList))); %numel(izpermislands.PixelIdxList)>0
        %                         % count sizes of clusters
        % %                         tempclustsizes = cellfun(@length,izpermislands.PixelIdxList);
        %                         % store size of biggest cluster
        %                         maxclustsize = max(cell2mat(cellfun(@numel,izpermislands.PixelIdxList)));
        % %                         [biggest,idx] = max(numPixels);
        % %                         BW(CC.PixelIdxList{idx})
        % %                         ixpc.MCmax_clust_size{ianimal}{introde} = [ixpc.MCmax_clust_size{ianimal}{introde}; max(tempclustsizes)];
        %                          ixpc.MCmax_clust_size{ianimal}{introde} = [ixpc.MCmax_clust_size{ianimal}{introde}; maxclustsize];
        %                     end
        %                 end
        
end
end