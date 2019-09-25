
% Similar to sj_getPFC_CA1_sigcorridxs
% Instead of glmfit file, use saved corr file, with data combined across epochs
% to get for each PFC cell, the corresponding CA1 cells with significant corrlns
% with sign

% ------
% From glm fit gather file, get for each PFC cell, the corresponding CA1 cells
% with significanct glm fit coreffs (or significant correlations) with sign

clear;
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
gatherdatafile = [savedir 'HP_ThetacovAndRipmodcorr_newtimefilter_alldata_gather3'];
load(gatherdatafile);


% Skip nan corrlns and no. of co-occurences <5/10 ?
rem = find(isnan(Sallr_epcomb));
rem2 = find(Sall_nsimul_epcomb<5);
allrem = union(rem, rem2);
Sallr_epcomb(allrem)=[]; Sallp_epcomb(allrem)=[]; runpairoutput_idx(allrem,:)=[]; 


% Now parse idxs to get sig cells for each PFC cell
% Idxs are [anim day CA1tet CA1cell PFCtet PFCcell]
% ----------------------------
% Initialize
cntcells = 0;
allPFCCA1sigidxs = [];
allPFCidxs = []; allCA1sigidxs = []; allCA1posidxs = []; allCA1negidxs = []; allsigr = []; allsigp = [];
% all unique indices of anim-day-PFCtet-PFCcell, collapsing across other idxs
uniqueIndices=unique(runpairoutput_idx(:,[1,2,5,6]),'rows');
% iterating only over the unique indices and finding matches in allindexes

for i=1:length(uniqueIndices)
    cntcells = cntcells+1;
    curridx = uniqueIndices(i,:);
    ind=find(ismember(runpairoutput_idx(:,[1 2 5 6]),curridx,'rows'))';
    
    % Corr data
    curr_rcorr = Sallr_epcomb(ind); curr_pcorr = Sallp_epcomb(ind);
    currCA1idxs = runpairoutput_idx(ind,[1 2 3 4]); % anim, day, CA1tet, CA1cell
    currsig = find(curr_pcorr<0.05);
    currsigpos = find(curr_rcorr(currsig)>0);
    currsigneg = find(curr_rcorr(currsig)<0);
    
    
    % Save for current PFC cell - both in struct, and variable format
    allPFCidxs(cntcells,:) = curridx;
    allCA1sigidxs{cntcells} = currCA1idxs(currsig,:);
    allCA1posidxs{cntcells} = currCA1idxs(currsig(currsigpos),:);
    allCA1negidxs{cntcells} = currCA1idxs(currsig(currsigneg),:);
    allsigr{cntcells} = curr_rcorr(currsig);
    allsigp{cntcells} = curr_pcorr(currsig);
    
    
    allPFCCA1sigcorridxs(cntcells).PFCidx = curridx;
    allPFCCA1sigcorridxs(cntcells).CA1sigidxs = allCA1sigidxs{cntcells};
    allPFCCA1sigcorridxs(cntcells).CA1posidxs = allCA1posidxs{cntcells};
    allPFCCA1sigcorridxs(cntcells).CA1negidxs = allCA1negidxs{cntcells};
    allPFCCA1sigcorridxs(cntcells).rsig = allsigr{cntcells};
    allPFCCA1sigcorridxs(cntcells).psig = allsigp{cntcells};
end

savefile = [savedir 'HP_allPFCCA1sigcorridxs'];
save(savefile,'allPFCCA1sigcorridxs');













