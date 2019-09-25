% From glm fit gather file, get for each PFC cell, the corresponding CA1 cells
% with significanct glm fit coreffs (or significant correlations) with sign

clear;
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
gatherdatafile = [savedir 'HP_ripmod_glmfit_gather2'];
load(gatherdatafile);

cnt=0; % How many kept? Any condition?
%Glm
allglmidxs = []; allmodelb=[]; allmodelp=[]; allmodelfits=[];
allnsig=[]; allnsigpos=[]; allnsigneg=[];  allfracsig=[]; allfracsigpos=[]; allfracsigneg=[];
%Corr
allcorridxs = []; allrcorr=[]; allpcorr=[]; allnsimul=[];

for an = 1:length(modf)
    for i=1:length(modf(an).output{1})
        % Check for empty output
        if ~isempty(modf(an).output{1}(i).nsig)
            cnt=cnt+1;
            anim_index{an}{cnt} = modf(an).output{1}(i).indices;
            % Only indexes
            %animindex=[an modf(an).output{1}(i).indices]; % Put animal index in front
            %allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
            %Sub indices for glm and corrcoef
            
            % Need to put "an" in front of idxs
            currglmidxs=[]; currcorridxs=[];
            ncurridxs = size(modf(an).output{1}(i).glmidxs,1); % No of idxs
            an_tmp = an*ones(ncurridxs,1);
            currglmidxs = [an_tmp, modf(an).output{1}(i).glmidxs];
            currcorridxs = [an_tmp, modf(an).output{1}(i).corridxs];
            
            allglmidxs = [allglmidxs; currglmidxs];
            allcorridxs = [allcorridxs; currcorridxs];
            % Glm data
            allmodelb = [allmodelb; modf(an).output{1}(i).allmodelb];
            allmodelp = [allmodelp; modf(an).output{1}(i).allmodelp];
            allmodelfits{cnt} = modf(an).output{1}(i).allmodelfits;
            modf(an).output{1}(i).nsig;
            allnsig = [allnsig; modf(an).output{1}(i).nsig'];
            allnsigpos = [allnsigpos; modf(an).output{1}(i).nsigpos'];
            allnsigneg = [allnsigneg; modf(an).output{1}(i).nsigneg'];
            allfracsig = [allfracsig; modf(an).output{1}(i).fracsig'];
            allfracsigpos = [allfracsigpos; modf(an).output{1}(i).fracsigpos'];
            allfracsigneg = [allfracsigneg; modf(an).output{1}(i).fracsigneg'];
            % Corr data
            allrcorr = [allrcorr; modf(an).output{1}(i).rcorr];
            allpcorr = [allpcorr; modf(an).output{1}(i).pcorr];
            allnsimul = [allnsimul; modf(an).output{1}(i).nsimul];
        end
    end
    
end

%keyboard;

% Skip bad fits, and maybe corrlns where no. of co-occurences are <10
rem = find( ((allmodelb>1) | (allmodelb<-1)) & allmodelp>0.99); % Corresponding p will usually be 1
%rem2 = find(allnsimul<10);
rem2=[];
allrem = union(rem, rem2);
allmodelb(allrem)=[]; allmodelp(allrem)=[]; allrcorr(allrem)=[]; allpcorr(allrem)=[]; allglmidxs(allrem,:)=[]; allcorridxs(allrem,:)=[];


% Now parse idxs to get sig cells for each PFC cell
% Idxs are [anim day epoch CA1tet CA1cell PFCtet PFCcell]
% ----------------------------
% Initialize
cntcells = 0;
allPFCCA1sigidxs = [];
allPFCidxs = []; allCA1sigidxs = []; allCA1posidxs = []; allCA1negidxs = []; allsigb = []; allsigp = [];
% all unique indices of anim-day-epoch-PFCtet-PFCcell, collapsing across other idxs
uniqueIndices=unique(allglmidxs(:,[1,2,3,6,7]),'rows');
% iterating only over the unique indices and finding matches in allindexes
for i=1:length(uniqueIndices)
    cntcells = cntcells+1;
    curridx = uniqueIndices(i,:);
    ind=find(ismember(allglmidxs(:,[1 2 3 6 7]),curridx,'rows'))';
    
    % Corr data
    currrcorr = allrcorr(ind); currpcorr = allpcorr(ind);
    % Glm data
    currp = allmodelp(ind); currb = allmodelb(ind);
    currCA1idxs = allglmidxs(ind,[1 2 3 4 5]); % anim, day, ep, CA1tet, CA1cell
    currsig = find(currp<0.05);
    currsigpos = find(currb(currsig)>0);
    currsigneg = find(currb(currsig)<0);
    
    
    % Save for current PFC cell - both in struct, and variable format
    allPFCidxs(cntcells,:) = curridx;
    allCA1sigidxs{cntcells} = currCA1idxs(currsig,:);
    allCA1posidxs{cntcells} = currCA1idxs(currsig(currsigpos),:);
    allCA1negidxs{cntcells} = currCA1idxs(currsig(currsigneg),:);
    allsigb{cntcells} = currb(currsig);
    allsigp{cntcells} = currp(currsig);
    
    
    allPFCCA1sigidxs(cntcells).PFCidx = curridx;
    allPFCCA1sigidxs(cntcells).CA1sigidxs = allCA1sigidxs{cntcells};
    allPFCCA1sigidxs(cntcells).CA1posidxs = allCA1posidxs{cntcells};
    allPFCCA1sigidxs(cntcells).CA1negidxs = allCA1negidxs{cntcells};
    allPFCCA1sigidxs(cntcells).bsig = allsigb{cntcells};
    allPFCCA1sigidxs(cntcells).psig = allsigp{cntcells};
end

savefile = [savedir 'HP_allPFCCA1sigidxs'];
save(savefile,'allPFCCA1sigidxs');

















