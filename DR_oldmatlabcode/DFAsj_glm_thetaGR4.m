function out = DFAsj_glm_thetaGR4(indices, excludetimes, spikes, cellinfo, varargin)
% Ver4 : Starting 10Feb2014 - Sync codes with everyone
% GLM model for theta times. Similar to DFAsj_glm_ripalign.m 
% Get the time series with theta times
% similar to the beginnin of the function"DFAsj_getthetacrosscov_timecondition.m" 


tetfilter = '';
% gideon
%excludetimes = [];
thrstime = 1; % Minimum length of includetime segments in sec
binsize = 0.2; % 200ms bins?

acrossregions = 0;

% lowsp_thrs = 5; %cm/sec
% highsp_thrs = lowsp_thrs;
% dospeed = 0;


for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'excludetimes'
            excludetimes = varargin{option+1};
        case 'thrstime'
            thrstime = varargin{option+1};  
        case 'binsize'
            binsize = varargin{option+1};
        case 'acrossregions'
            acrossregions = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

% First, get cell indices for CA1 and PFC
% ----------------------------------------

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
    

% Get epochtimes, and implement the timecondition
% ------------------------------------------------

ind = cellsi(1,:); % day.ep.tet.cell for 1st CA1 cell 
ttemp = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);

totaleptime = diff(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange)./10000; % in secs
excltime = sum(diff(excludetimes'));
% Use Include periods
% --------------------------
% Get IncludeTimes from flanking edges in excludetimes and epoch start-end
epstend = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange./10000;
incl=[];
incl(:,1) = excludetimes(1:end-1,2);
incl(:,2) = excludetimes(2:end,1);
incl = [epstend(1),incl(1,1) ;incl];
incl = [incl; incl(end,2),epstend(2)];

% Length of include periods
incl_lths = diff(incl')';
% Discard anything < thrstime
discard = find(incl_lths<thrstime);
incl(discard,:)=[];


% Now, "incl" has a set of start and endtimes for "velid" timeperiods which are atleast "thrstime" long
% For each cell's spikes, the spiketimes have to be runthrough isIncluded as: 
% t1inc = t1(find(isIncluded(t1, incl)));
% where t1 is the spiketimes for the cell in the epoch. eg. t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
% ------------------------------------------------------------------------



% ------------------
% Parse spikes 
for i=1:size(cellsi,1)
    i;
    eval(['spikesi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
    % gideon added, verify with shantanu
    eval(['spikesiinc{',num2str(i),'}= spikesi{',num2str(i),'}(find(isIncluded(spikesi{',num2str(i),'},incl)))'])
end

for i=1:size(cellsp,1)
    i;
    eval(['spikesp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
 % gideon added, verify with shantanu
    eval(['spikespinc{',num2str(i),'}= spikesp{',num2str(i),'}(find(isIncluded(spikesp{',num2str(i),'},incl)))'])

end

% Going over all incl intervals, and from each interval extracting bin
% times of 200ms going from the interval start, ending >=200ms from its end
% thetaBins will hold the start times of all 200ms bins
thetaBins=[];
thetaBinSize=0.2;
for ii=1:length(incl),
    curIntervalSize=incl(ii,2)-incl(ii,1);
    numCurBins=floor(curIntervalSize/thetaBinSize);
    curBins=incl(ii,1):0.2:(incl(ii,1)+0.2*(numCurBins-1));
    thetaBins=[thetaBins curBins];
end




% ---------------------------------- ----------------------------------
% HAS TO BE UPDATED FROM THIS POINT ON
% DONE
% ---------------------------------- ----------------------------------

spikeBinsHC=[];
spikeBinsPFC=[];
for k=1:nCA1cells
    
spikesToBins=lookup(spikesiinc{k},thetaBins,-1);
tabulateSpikes=tabulate(spikesToBins);
spikesInBins=NaN(1,length(thetaBins));
%if ~isempty(tabulateSpikes)
spikesInBins(1:length(tabulateSpikes))=tabulateSpikes(:,2);
%end
spikeBinsHC=[spikeBinsHC spikesInBins'];
end

for k=1:nPFCcells
    
spikesToBins=lookup(spikespinc{k},thetaBins,-1);
tabulateSpikes=tabulate(spikesToBins);
spikesInBins=NaN(1,length(thetaBins));
%if ~isempty(tabulateSpikes)
spikesInBins(1:length(tabulateSpikes))=tabulateSpikes(:,2);
%end
spikeBinsPFC=[spikeBinsPFC spikesInBins'];
end


Xmat=spikeBinsHC;

% Get glmfits for the model using Poisson distr, which uses log link function by default:
% PFC trialresps(y), X = constant and CA1 trialResps, (log(u)) = Xb



% No need to add a column of constant at beginning. glm automatically does it
% nobs = size(Xmat,1); constants = ones(nobs,1); Xmat = [constants, Xmat];

allmodelfits = []; allmodelb = []; allmodelp=[]; allidxs = []; allidxs2=NaN(20,100);
allmodelb2=NaN(20,100);
nsig = []; nsigpos = []; nsigneg = []; fracsig =[]; fracsigpos=[]; fracsigneg=[];
% For each PFC cell
for i=1:size(cellsp,1)
    cellsp(i,:);
   % eval(['y = ripplemodp{',num2str(i),'};']);
   y=spikeBinsPFC(:,i); 
   [b, ~, stats] = glmfit(Xmat,y,'poisson');
    allmodelfits(i).stats = stats;
    allmodelb = [allmodelb; b(2:end)]; % Coefficients. Same as stats.beta
    allmodelp = [allmodelp; stats.p(2:end)];  
    
    % Save index as [day epoch CA1tet CA1cell PFCtet PFCcell]
    PFCdet = repmat([cellsp(i,3), cellsp(i,4)], nCA1cells,1 );
    curridxs = [cellsi,PFCdet]; 
    allidxs = [allidxs; curridxs];   
    
    % Indices in new form. Each line is one ensemble of cells in the
    % following form:
    % day,epoch,hctet,hccell,hctet,hccell...,hctet,hccell,pfctet,pfccell
    hctetcell=curridxs(:,3:4);
    hctetcell2=reshape(hctetcell',1,length(hctetcell)*2);
    curridxs2=[curridxs(1,1:2) hctetcell2 curridxs(1,end-1:end)];
    allidxs2(i,1:length(curridxs2))=curridxs2;
    allmodelb2(i,1:length(b))=b;
    % For each PFC neuron, what fraction of CA1 cells were significantly predictive - positive or negative?
    currsig = find(stats.p(2:end) < 0.05);
    nsig(i) = length(currsig);
    fracsig(i) = nsig(i)./nCA1cells;
    bsig = b(currsig+1);
    nsigpos(i) = length(find(bsig>0)); fracsigpos(i) = nsigpos(i)./nCA1cells;
    nsigneg(i) = length(find(bsig<0)); fracsigneg(i) = nsigneg(i)./nCA1cells;
    
end
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
out.glmidxs = allidxs;
out.glmidxs2=allidxs2;
out.allmodelb2=allmodelb2;
out.allmodelb = allmodelb;
out.allmodelp = allmodelp;
out.allmodelfits = allmodelfits; % has allmodelfits.stats, which has all infor about model
%try
    out.nsig = nsig; % No. of sig CA1 cells in model for each corresponding PFC cell.
% catch
%     day, epoch
%     keyboard;
% end
out.nsigpos = nsigpos;
out.nsigneg = nsigneg;
out.fracsig = fracsig;
out.fracsigpos = fracsigpos;
out.fracsigneg = fracsigneg;

% Corr coeff
%out.corridxs = corridxs;
%out.rcorr = rcorr;
%out.pcorr = pcorr;
%out.nsimul = nsimul; % No. of co-occurences for each pair





