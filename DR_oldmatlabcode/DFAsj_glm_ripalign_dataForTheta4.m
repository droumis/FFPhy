function out = DFAsj_glm_ripalign_dataForTheta4(indices, excludetimes, ripplemod, cellinfo, varargin)
% Ver4 : Starting 10Feb2014 - Sync codes with everyone
% out = DFAsj_getripalignspiking(spike_index, excludeperiods, spikes, ripples, tetinfo, options)

% See also sj_glm_ripalign1

% Need ripplemod only. IF across regions, need cellinfo
% Need last spikes, ripples, pos, linpos, only if calculating ripplealign from scratch as follows
% See DFSsj_getripalignspiking
% Use tetinfo and tetfilter passed in, or redefine here to get riptets
% Then use ripples to getriptimes. Use inter-ripple-interval of 1 sec, and use a low-speed criterion.
% Then align spikes to ripples


tetfilter = '';
excludetimes = [];
maxcell = 0;
minstd = 3;
acrossregions = 0;
% lowsp_thrs = 5; %cm/sec
% highsp_thrs = lowsp_thrs;
% dospeed = 0;
% 
% % For ripple trigger
% % ------------------
% binsize = 10; % ms
% pret=550; postt=550; %% Times to plot
% push = 500; % either bwin(2) or postt-trim=500. For jittered trigger in background window
% trim = 50;
% cellcountthresh = 3;  % Can be used to parse ripples
% smwin=10; %Smoothing Window - along y-axis for matrix. Carry over from ...getrip4
% 
% rwin = [0 150];
% bwin = [-500 -300];
% push = 500; % either bwin(2) or postt-trim=500. If doing random events. See ... getrip4


for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'excludetimes'
            excludetimes = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'dospeed'
            dospeed = varargin{option+1};
        case 'acrossregions'
            dospeed = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

day = indices(1,1);
epoch = indices(1,2);

cellsi=[]; usecellsi=0; % CA1 cells
cellsp=[]; usecellsp=0; % PFC cells

totalcells = size(indices,1);
for i=1:totalcells
    currarea = cellinfo{indices(i,1)}{indices(i,2)}{indices(i,3)}{indices(i,4)}.area;
    if strcmp(currarea,'PFC'),
        cellsp = [cellsp; day epoch indices(i,3) indices(i,4)];
        usecellsp = usecellsp+1;
    else
        cellsi = [cellsi; day epoch indices(i,3) indices(i,4)];
        usecellsi = usecellsi+1;
    end
end
nCA1cells = size(cellsi,1); nPFCcells = size(cellsp,1);
    
% ------------------
% Parse ripplemoddata - don;t really need to do this, but following format from sj_glm_ripalign1
for i=1:size(cellsi,1)
    i;
    eval(['ripplemodi{',num2str(i),'}= ripplemod{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.trialResps;']);
end

for i=1:size(cellsp,1)
    i;
    eval(['ripplemodp{',num2str(i),'}= ripplemod{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.trialResps;']);
end  
    

% Get glmfits for the model using Poisson distr, which uses log link function by default:
% PFC trialresps(y), X = constant and CA1 trialResps, (log(u)) = Xb

% X will be the same for all PFC cells
Xmat = [];
for i=1:size(cellsi,1)
    eval(['currResps = ripplemodi{',num2str(i),'};']);
    Xmat = [Xmat, currResps]; % Rows are observations, Columns are triaslResps for current CA1 cell
end
Ymat=[];
for i=1:size(cellsp,1)
     eval(['y = ripplemodp{',num2str(i),'};']);
     Ymat=[Ymat y];
end
% No need to add a column of constant at beginning. glm automatically does it
% nobs = size(Xmat,1); constants = ones(nobs,1); Xmat = [constants, Xmat];

allidxs2=NaN(20,100);
% For each PFC cell
for i=1:size(cellsp,1)
%     cellsp(i,:);
%     eval(['y = ripplemodp{',num2str(i),'};']);
%     [b, ~, stats] = glmfit(Xmat,y,'poisson');
%     allmodelfits(i).stats = stats;
%     allmodelb = [allmodelb; b(2:end)]; % Coefficients. Same as stats.beta
%     allmodelp = [allmodelp; stats.p(2:end)];  
    
    % Save index as [day epoch CA1tet CA1cell PFCtet PFCcell]
    PFCdet = repmat([cellsp(i,3), cellsp(i,4)], nCA1cells,1 );
    curridxs = [cellsi,PFCdet];
%     allidxs = [allidxs; curridxs];   
    % Indices in new form. Each line is one ensemble of cells in the
    % following form:
    % day,epoch,hctet,hccell,hctet,hccell...,hctet,hccell,pfctet,pfccell
    hctetcell=curridxs(:,3:4);
    hctetcell2=reshape(hctetcell',1,length(hctetcell)*2);
    curridxs2=[curridxs(1,1:2) hctetcell2 curridxs(1,end-1:end)];
    allidxs2(i,1:length(curridxs2))=curridxs2;
end
    % For each PFC neuron, what fraction of CA1 cells were significantly predictive - positive or negative?
%     currsig = find(stats.p(2:end) < 0.05);
%     nsig(i) = length(currsig);
%     fracsig(i) = nsig(i)./nCA1cells;
%     bsig = b(currsig+1);
%     nsigpos(i) = length(find(bsig>0)); fracsigpos(i) = nsigpos(i)./nCA1cells;
%     nsigneg(i) = length(find(bsig<0)); fracsigneg(i) = nsigneg(i)./nCA1cells;
%     
%     numRips=size(Xmat,1);
%     allErrReal=[];
%     allErrShuf=[];
%     allFluctErr=[];
%     allFluctErrShuf=[];
%     for ii=1:100
%         ripidxs=randperm(numRips);
%         dataPercentForTrain=0.9;
%         Ntrain=ripidxs(1:round(numRips*dataPercentForTrain));
%         Ntest=ripidxs(round(numRips*dataPercentForTrain)+1:numRips);
%         [btr, ~, statstr] = glmfit(Xmat(Ntrain,:),y(Ntrain),'poisson');
%         % btr(statstr.p > 0.05)=0;
%         yfit = glmval(btr, Xmat(Ntest,:),'log',statstr,0.95);
%         Ntestshufd=Ntest(randperm(length(Ntest)));
%         %  yfitshuf = glmval(btr, Xmat(Ntestshufd,:),'log',statstr,0.95);
%         errReal=mean((yfit-y(Ntest)).^2);
%         %errShuf=mean((yfitshuf-y(Ntest)).^2);
%         errShuf=mean((yfit-y(Ntestshufd)).^2);
%         
%         allErrReal=[allErrReal errReal];
%         allErrShuf=[allErrShuf errShuf];
%         fluctErr=mean(abs((yfit>mean(yfit))-(y(Ntest)>mean(y(Ntest)))));
%         fluctErrShuf=mean(abs((yfit>mean(yfit))-(y(Ntestshufd)>mean(y(Ntestshufd)))));
%         allFluctErr=[allFluctErr fluctErr];
%         allFluctErrShuf=[allFluctErrShuf fluctErrShuf];
%         
%     end
%     allErrRealAllCells=[allErrRealAllCells mean(allErrReal)];
%     allErrShufAllCells=[allErrShufAllCells mean(allErrShuf)];
%     [r1,p1]=ttest(allErrReal,allErrShuf);
%     allSigPrediction=[allSigPrediction p1];
%     allFluctErrAllCells=[allFluctErrAllCells mean(allFluctErr)];
%     allFluctErrShufAllCells=[allFluctErrShufAllCells mean(allFluctErrShuf)];
%     
%     %     [a1 b1]=hist(allErrReal,0:0.1:3);
%     %     [a2 b2]=hist(allErrShuf,0:0.1:3);
% %     figure;plot(b1,a1);hold on;plot(b2,a2,'r')
% %     legend('real data errors','shuffled data errors')
% %     [r1 p1]=ttest(allErrReal,allErrShuf)
% %     errorDecrease=mean(allErrReal)/mean(allErrShuf);
% % title(['ttest p-val=' num2str(p1) ' , error decrease= ' num2str(errorDecrease)]); 
% %     
%     %keyboard
% end
% 
% % Get pairwise correlations between CA1 and PFC - to compare to model coeffs
% % Skip getting significance from shuffling for now.
% corridxs = []; rcorr = []; pcorr = []; nsimul = [];
% for pp=1:size(cellsp,1)
%      eval(['y = ripplemodp{',num2str(pp),'};']); % PFC cell
%      for ii=1:size(cellsi,1)
%          eval(['x = ripplemodi{',num2str(ii),'};']); % CA1 cell
%          [r, p] = corrcoef(x,y);
%          rcorr = [rcorr; r(1,2)];
%          pcorr = [pcorr; p(1,2)];
%          corridxs = [corridxs; day epoch cellsi(ii,3) cellsi(ii,4) cellsp(pp,3) cellsp(pp,4)];
%          
%          % Get number of "trials/ripples" with co-occurences as well
%          nsimul = [nsimul; length(find((x~=0) & (y~=0)))];
%      end
% end

% Leave this for the calling script
% Skip bad fits, and corrlns. where no. of co-occurences are <10
% rem = find( (allmodelb>1) | (allmodelb<-1)); % Corresponding p will be 1
% rem2 = find(nsimul<10)
% allrem = union(rem, rem2);
% allmodelb(allrem)=[]; allmodelp(allrem)=[]; rcorr(allrem)=[]; pcorr(allrem)=[]; allidxs(allrem,:)=[]; corridxs(allrem,:)=[];
% sigglm = find(allmodelp < 0.05);
% sigcorr = find(pcorr < 0.05);


% % ------ 
% % Output
% % ------
out.indices = indices;
% Glm fit
%out.glmidxs = allidxs;
out.glmidxs2=allidxs2;
out.Xmat=Xmat;
out.Ymat=Ymat;
% 
% out.allmodelb = allmodelb;
% out.allmodelp = allmodelp;
% out.allmodelfits = allmodelfits; % has allmodelfits.stats, which has all infor about model
% %try
%     out.nsig = nsig; % No. of sig CA1 cells in model for each corresponding PFC cell.
% % catch
% %     day, epoch
% %     keyboard;
% % end
% out.nsigpos = nsigpos;
% out.nsigneg = nsigneg;
% out.fracsig = fracsig;
% out.fracsigpos = fracsigpos;
% out.fracsigneg = fracsigneg;
% % Corr coeff
% out.corridxs = corridxs;
% out.rcorr = rcorr;
% out.pcorr = pcorr;
% out.nsimul = nsimul; % No. of co-occurences for each pair
% out.allErrReal=allErrRealAllCells;
% out.allErrShuf=allErrShufAllCells;
% out.predictSig=allSigPrediction;
% 
% out.allFluctErrReal=allFluctErrAllCells;
% out.allFluctErrShuf=allFluctErrShufAllCells;
% 
% 
