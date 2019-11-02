%load('bond04_vecLF_fold.mat');
%load('bond04_vecLF.mat');

day = 4;
ep = 2;

animalname = 'Bond';
    animalinfo = animaldef(animalname);
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',day);

eventoption = 2;  %  == 2 for ripplescons, == 3 for WG
extratime = 500;  % ms of extra time before and after ripple to decode 
WGexclusion = 1;
if eventoption == 2
    cons_name1 = 'ripplescons';
    tetfilter1 = 1;     % 1 for validripples CA1, 2 for CA3-DG
    consensus_numtets_ripc = 3;
    minthresh_ripc = 7;
    exclusion_dur_ripc = 0;
    minvelocity_ripc = 0;
    maxvelocity_ripc = 4;  
elseif eventoption == 3
    cons_name1 = 'gammafcons';
    tetfilter1 = 2;     % 1 for validripples CA1, 2 for CA3-DG
    consensus_numtets_ripc = 2;
    minthresh_ripc = 2;
    exclusion_dur_ripc = 0;
    minvelocity_ripc = 8;
    maxvelocity_ripc = inf;      
end



if WGexclusion
    cons_name2 = 'gammafcons';
    tetfilter2 = 2;     % 1 for validripples CA1, 2 for CA3-DG
    consensus_numtets_gfc = 2;
    minthresh_gfc = 2;
    exclusion_dur_gfc = 0;
    minvelocity_gfc = 0;
    maxvelocity_gfc = inf;
        exclusion2_name = 'ripplescons';
        exclusion2_tetfilter = 1;
        exclusion2_consensus_numtets = 3;
        exclusion2_minthresh = 2;
        exclusion2_exclusion_dur = 0;
        exclusion2_maxvelocity = 4;
    
    output2 = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],cons_name2,tetfilter2,...
        'consensus_numtets',consensus_numtets_gfc,'minthresh',minthresh_gfc,...
        'exclusion_dur',exclusion_dur_gfc,'minvelocity',minvelocity_gfc,'maxvelocity',maxvelocity_gfc,...
        'exclusion2',exclusion2_name,exclusion2_tetfilter,exclusion2_consensus_numtets,exclusion2_minthresh,...
            exclusion2_exclusion_dur,exclusion2_maxvelocity)
    consvec2 = output2{day}{ep}.cons;
    consvectimes2 = output2{day}{ep}.time;
    gfperiods = vec2list(consvec2,consvectimes2);
        gfdur = sum(gfperiods(:,2)-gfperiods(:,1));
        disp(sprintf('%d s of wave gamma detected',round(gfdur)))
    
end

load('/opt/data13/kkay/Bon/bonlinpos04.mat');


%%% basic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
postimevec=linpos{4}{2}.statematrix.time;            
poslin=vecLF(:,2);  
    poslin2 = poslin';

%%% basic values:   xdel, xbins, dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdel=0.1; 
xbins=min(poslin):xdel:max(poslin);        % -3 : 0.1 : 3
    numlinbins = length(xbins);
dt=postimevec(2)-postimevec(1);         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make State M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stateM=[];
[~,binind]=histc(poslin,xbins);        % histogram of linearized positions
xbinstate=[binind(1:end-1) binind(2:end)];  % [ <current step lin pos>  <next step lin pos>  ]
%by column = departuring
for bb = 1:numlinbins
    % nextstates is the next states;
    nextstates = xbinstate( xbinstate(:,1)==bb , 2  ); 
    if ~isempty(nextstates)
        % histograms the nextstates and normalizes the distribution (why?)
        stateM(:,bb) = histc( nextstates , linspace(1,numlinbins,numlinbins))./size(nextstates,1 );
    elseif isempty(nextstates)
        stateM(:,bb) = zeros(1,numlinbins);
    end
end
K = inline('exp(-(x.^2+y.^2)/2/sig^2)');                        %gaussian
[dx,dy] = meshgrid([-1:1]);
sig = 0.5;
weight = K(sig,dx,dy)/sum(sum(K(sig,dx,dy)));                   %normalizing weights
stateM_gaus = conv2(stateM,weight,'same');                      %gaussian smoothed
stateM_gausnorm = stateM_gaus*diag(1./sum(stateM_gaus,1));      %normalized to confine probability to 1




%%%%%%  encode per tetrode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mark0 = {};       % mark0_t1  
spktimes = {};    % time_t1
spktimes2 = {};   % time2_t1
procInd0 = {};    % procInd0   % pos time raster of spikes ??
procInd1 = {};    % procInd1   % pos ind assignment of spikes
procInd = {};
spikeT = {};
spike = {};
rawInd0 = {};
markAll = {};
Xnum = {};
Lint = {};


for tet = selected_tets
    
    daydir = getdaydir(animalname);
    filedata = loadparamsfile(daydir,day,tet);
    
    
    % identify spikes to use in encoding
    
    % spikes in experimenter-transcribed epoch
    inds1 =  (  filedata.params(:,1)/10000  >=  postimevec(1)  )  &  ( filedata.params(:,1)/10000 <= postimevec(end) );
    % spikes that has 100 uV amplitude in at least one channel 
    inds2 = (filedata.params(:,2) > 100) | (filedata.params(:,3) > 100) | (filedata.params(:,4) > 100) | (filedata.params(:,5) > 100);
    % spikes that occur above 6 cm/s
    spikeposinds = lookup(filedata.params(:,1)/10000,pos{day}{ep}.data(:,1));
    spikevels = pos{day}{ep}.data(spikeposinds,5);
    inds3 = spikevels > 6;
    
    
    % determine final set of inds
    if ~WGexclusion
        inds = inds1 & inds2 & inds3;   
    else
        % spikes that occur during WG periods
        inds4 = isExcluded(filedata.params(:,1)/10000,gfperiods);
        % intersection: valid spike inds
        inds = inds1 ;%& inds2 & inds3 & ~inds4;  % spikes to encode with
        testinds = inds1;% & inds2 & inds3;   % without WG filter
        indsALL = inds1; % & inds2;  % inclusive, to decode with
        disp(sprintf('%d of %d spikes excluded for being in WG period',sum(testinds)-sum(inds),sum(testinds)))          
    end

    % mark vector (4 channel amplitudes)
    mark0{tet}=[filedata.params(inds,2) filedata.params(inds,3) filedata.params(inds,4)  filedata.params(inds,5)];
    mark0_ALL{tet} = [filedata.params(indsALL,2) filedata.params(indsALL,3) filedata.params(indsALL,4)  filedata.params(indsALL,5)];
    
    % spike timestamps  (inclusive, to decode)
    spktimes{tet}=filedata.params(indsALL,1);   
    
    % spike timestamps in sec  (to encode)
    spktimes2{tet}=filedata.params(inds,1)/10000;
    
    % procInd0
    [~,procInd1_ALL{tet}]=histc(spktimes{tet}/10000,postimevec);
    
    [procInd0{tet},procInd1{tet}]=histc(spktimes2{tet},postimevec);
    procInd{tet}=find(procInd0{tet});
    spikeT{tet}=postimevec(procInd{tet});
    spike{tet}=procInd0{tet}';
    [~,rawInd0{tet}]=histc(spktimes2{tet},spktimes2{tet});
    markAll{tet}(:,1)=procInd1{tet};
    markAll{tet}(:,2:5)=mark0{tet}(rawInd0{tet}(rawInd0{tet}~=0),:);
    mdel=20; ms = min(min(markAll{tet}(:,2:5))):mdel:max(max(markAll{tet}(:,2:5)));
    sxker=2*xdel; smker = mdel; T=size(postimevec,1);
    occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
    %occ: columns are identical; occupancy based on position; denominator
    
    % Xnum: For each spike (each columns), calculates a position gaussian
    % centered at position at which spike occurred
    Xnum{tet}=normpdf(xbins'*ones(1,length(spktimes2{tet})),...
                      ones(length(xbins),1)*poslin2(procInd1{tet}),...   % column of identical MUs (linear position at which spike occurred)
                      sxker);
    %Xnum: Gaussian kernel estimators for position
    Lint{tet} = sum(Xnum{tet},2)./occ(:,1)./dt; %integral
    Lint{tet} = Lint{tet}./sum(Lint{tet});
    
end




% Collect all spikes, marks, and procInd, and lengths
spktimes_all = [];          % formerly time0 = [];   % all spike times
spktimes_all_sorted = [];   % formerly time = [];    % sorted spike times
sinds = [];                 % formerly timeInd = [];  % all sorted spike inds
mark0_all = [];
procInd1_all = [];
len = {};
for tet = selected_tets
    
    numspktet = length(spktimes{tet});
    spktimes_all = [spktimes_all ;  spktimes{tet}    tet * ones(numspktet,1) ];   % * tag with the tet number
    mark0_all=[ mark0_all ; mark0_ALL{tet} ];
    procInd1_all=[procInd1_all ; procInd1_ALL{tet}];
    len{tet} = length(spktimes{tet});
    
end

% Sort all spikes and marks
[spktimes_all_sorted,sinds] = sort(spktimes_all(:,1));
mark0_all=mark0_all(sinds,:);
procInd1_all=procInd1_all(sinds,:);
 
% Construct tet_ind "indicator matrix" :    [  <time point>  x  <which tetrode spikes> ]

tet_ind = zeros(length(spktimes_all_sorted),max(selected_tets));

for i=1:length(spktimes_all)
    parenttet = spktimes_all(sinds(i),2);   % tet responsible for the spike 
    tet_ind(i,parenttet)=1; 
end

tet_sum=tet_ind.*cumsum(tet_ind,1); %row: time point; column: index of spike per tetrode


%%
if 1
    mdel=10; ms=min(mark0_all(:)):mdel:max(mark0_all(:));
else
   mdel=20; ms=100:mdel:max(mark0_all(:));  
end
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_all=normpdf(xbins'*ones(1,length(spktimes_all)),ones(length(xbins),1)*poslin2(procInd1_all),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_all=sum(Xnum_all,2)./occ(:,1)./dt; %integral
Lint_all=Lint_all./sum(Lint_all);
%Lint: conditional intensity function for the unmarked case



%%%% Event identification (ripples or WG) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if eventoption == 1
    
    load('/opt/data13/kkay/Bon/bonripples04.mat');
    %selected_tets=[1 2 4 5 7 10 11 12 13 14 17 18 19 20 22 23 27 29];
    ex=4;ep=2;
    
    ti=round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector
    %position time; from 2461011 to 3404979; length(ti)~~max(spike_tim)
    
    % Make ind_ripple 
    ind_ripple=zeros(size(ti,2),size(selected_tets,2));  % <time bin> x <tetnum>
    % iterate tetrodes
    ripple_all = [];
    for ii=1:size(selected_tets,2)
        k=selected_tets(ii);  % tet
        ripple_s=round(ripples{ex}{ep}{k}.starttime*1000);
        ripple_e=round(ripples{ex}{ep}{k}.endtime*1000);
        for l=1:size(ripple_s,1)
            ripple=[ripple_s(l) ripple_e(l)];
            ripple_all=[ripple_all;ripple];
        end
        clear ripple_s ripple_e ripple;
        for s=1:size(ripple_all,1)
            ind_ripple( find(  ti>=ripple_all(s,1)   &   ti<=ripple_all(s,2)  ),ii )=1;
        end
        ripple_all=[];
    end
    
    ripple_acr=zeros(size(ind_ripple,1),1);
    %ripple_acr=sum(ind_ripple,2);
    ripple_acr(   find(  sum(ind_ripple,2)  )   )=1;
    %%%indicator for ripples pooled across tetrodes
    
    %ripple_acr(find(sum(ind_ripple,2)>=6))=1;
    %%%indicator for ripples pooled across tetrodes where more than 6/9
    %%%tetrodes have spiking   ************** erroneous
    
    clear startT endT
    % iterate through 1-ms bin times
    for k=1:size(ti,2)-1
        if ripple_acr(k)==0 && ripple_acr(k+1)==1
            startT(k+1)=k+1;
        elseif ripple_acr(k)==1 && ripple_acr(k+1)==0
            endT(k)=k;
        end
    end
    if ripple_acr(size(ti,2))==1
        endT(size(ti,2))=size(ti,2);
    end
    
    ripple_seg=[find(startT)' find(endT)']; %index for ripple segments


elseif eventoption == 2
    
    output1 = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],cons_name1,tetfilter1,...
        'consensus_numtets',consensus_numtets_ripc,'minthresh',minthresh_ripc,...
        'exclusion_dur',exclusion_dur_ripc,'minvelocity',minvelocity_ripc,'maxvelocity',maxvelocity_ripc);
    consvec = output1{day}{ep}.cons;
    consvectimes = output1{day}{ep}.time;
    riplist = vec2list(consvec,consvectimes);
    
    disp(sprintf('%d Ripples detected',size(riplist,1)))
    
    ti=round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector (from position time)
    
    ind_ripple=list2vec(riplist,ti/1000);  % <time bin> x <consvec>
    
    ripple_acr=zeros(size(ind_ripple,1),1);
    %ripple_acr=sum(ind_ripple,2);
    ripple_acr(   find(  sum(ind_ripple,2)  )   ) = 1;
    %%%indicator for ripples pooled across tetrodes
    
    %ripple_acr(find(sum(ind_ripple,2)>=6))=1;
    %%%indicator for ripples pooled across tetrodes where more than 6/9
    %%%tetrodes have spiking   ************** erroneous
    
    clear startT endT
    % iterate through 1-ms bin times
    for k=1:size(ti,2)-1
        if ripple_acr(k)==0 && ripple_acr(k+1)==1
            startT(k+1)=k+1;
        elseif ripple_acr(k)==1 && ripple_acr(k+1)==0
            endT(k)=k;
        end
    end
    if ripple_acr(size(ti,2))==1
        endT(size(ti,2))=size(ti,2);
    end
    
    ripple_seg=[find(startT)' find(endT)']; %index for ripple segments    
    
elseif eventoption == 3
    
        output1 = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],cons_name1,tetfilter1,...
        'consensus_numtets',consensus_numtets_ripc,'minthresh',minthresh_ripc,...
        'exclusion_dur',exclusion_dur_ripc,'minvelocity',minvelocity_ripc,'maxvelocity',maxvelocity_ripc);
    consvec = output1{day}{ep}.cons;
    consvectimes = output1{day}{ep}.time;
    riplist = vec2list(consvec,consvectimes);
    
    disp(sprintf('%d Wave Gamma episodes detected',size(riplist,1)))
    
    ti=round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector (from position time)
    
    ind_ripple=list2vec(riplist,ti/1000);  % <time bin> x <consvec>
    
    ripple_acr=zeros(size(ind_ripple,1),1);
    %ripple_acr=sum(ind_ripple,2);
    ripple_acr(   find(  sum(ind_ripple,2)  )   ) = 1;
    %%%indicator for ripples pooled across tetrodes
    
    %ripple_acr(find(sum(ind_ripple,2)>=6))=1;
    %%%indicator for ripples pooled across tetrodes where more than 6/9
    %%%tetrodes have spiking   ************** erroneous
    
    clear startT endT
    % iterate through 1-ms bin times
    for k=1:size(ti,2)-1
        if ripple_acr(k)==0 && ripple_acr(k+1)==1
            startT(k+1)=k+1;
        elseif ripple_acr(k)==1 && ripple_acr(k+1)==0
            endT(k)=k;
        end
    end
    if ripple_acr(size(ti,2))==1
        endT(size(ti,2))=size(ti,2);
    end
    
    ripple_seg=[find(startT)' find(endT)']; %index for ripple segments    
    
end


%%
ripple_dur=ripple_seg(:,2)-ripple_seg(:,1);
%1203, 1237

xi=round(time/10);  %spike times
% tot_ripple=zeros(length(ripple_dur),1);
% for i=1:length(ripple_dur)
%     numSteps=ripple_dur(i);
%     spike_tim=ripple_seg(i,1):ripple_seg(i,2); 
%     a1=0;
%     for t=1:numSteps
%         tt=spike_tim(t);
%         a=find(xi==ti(tt));
%         a1=a1+length(a);
%     end
%     tot_ripple(i)=a1;
% end


numrips = size(ripple_seg,1);

for ripnum = 1:numrips

    % ms time bins to decode
    spike_tim = (ripple_seg(ripnum,1)  - extratime ) : (ripple_seg(ripnum,2) + extratime );        %from 1 to 90000~

    % ms time of ripple start
    ripms_start = lookup(ripple_seg(ripnum,1),spike_tim);
    ripms_end = lookup(ripple_seg(ripnum,2),spike_tim);
    
    % clock times for ripple
    rip_starttime = postimevec(1) + ripple_seg(ripnum,1)/1000;
    rip_endtime = postimevec(1) + ripple_seg(ripnum,2)/1000;
    
    % position of the animal during ripple (should barely change)
    PADDING = 20;  % how many position samples (29.97 Hz) before and after to collect pos for
    rippos = [];  % vector of positions over the course of the ripple
        rippos_startind = lookup(rip_starttime,postimevec) - PADDING;
        rippos_endind = lookup(rip_endtime,postimevec) + PADDING;
        
    if rippos_startind <= 0 || rippos_endind > length(postimevec)
        disp('cant plot event at edge of epoch')
        continue
    end
        
    rippos = poslin(rippos_startind:rippos_endind);
    rippos_timevec = 1000 * ( postimevec(rippos_startind:rippos_endind) - rip_starttime );
    
    
%%
clear onstep postx postxM_r;
n=size(stateM_gausnorm,1);              % # of position; length(xbins)
postx=ones(n,1)./n;                     % uniform prior
numSteps=length(spike_tim);
postxM_r=zeros(n,numSteps);

spike_r=zeros(max(selected_tets),numSteps);
xi=round(spktimes_all_sorted/10);       %spike times

tic
for t=1:numSteps
    tt=spike_tim(t);
    onestep=stateM_gausnorm*postx;
    L=ones(size(postx));
    a=find(xi==ti(tt));
    
    if isempty(a)==1 %if no spike occurs at time t
        L=exp(-Lint_all.*dt);
    elseif isempty(a)==0 %if spikes
        l=zeros(n,length(a));
        for j=1:length(a)
            jj=a(j);
            tet=find(tet_ind(jj,:));
            spike_r(tet,t)=1;
            i=tet_sum(jj,tet);
            l0=normpdf(markAll{tet}(i,2)*ones(1,length(spktimes2{tet})),markAll{tet}(:,2)',smker)...
                .*normpdf(markAll{tet}(i,3)*ones(1,length(spktimes2{tet})),markAll{tet}(:,3)',smker)...
                .*normpdf(markAll{tet}(i,4)*ones(1,length(spktimes2{tet})),markAll{tet}(:,4)',smker)...
                .*normpdf(markAll{tet}(i,5)*ones(1,length(spktimes2{tet})),markAll{tet}(:,5)',smker);
            l1=Xnum{tet}*l0'./occ(:,1)./dt;
            l2=l1.*dt.*exp(-Lint{tet}.*dt);
            l2=l2./sum(l2);
            l(:,j)=l2;
        end
        L=prod(l,2);
        L=L./sum(L);
    end
    postx = onestep.*L./sum(onestep.*L);
    postxM_r(:,t)=postx;
    clear onestep a l L;
end
toc

spike_r_raw = spike_r;



%%




if 0
    
    figure;
    subplot(3,2,[1 3 5]);
    imagesc(spike_r);colormap(flipud(gray));freezeColors
    set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
    xlabel('Time (ms)');ylabel('Tetrode');
    title(['ripple = ',num2str(ripnum)]);
    box off
    
    ind_seg1=1:21;
    ind_seg3=41:61;
    
    subplot(3,2,2);
    imagesc(1:numSteps,xbins(ind_seg1),postxM_r(ind_seg1,:));
    title('postx, clusterless','FontSize',18);
    ylabel('Linearized position','FontSize',14);
    set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
    colormap(flipud(hot(256)));
    caxis([0 0.3]);
    xlim([0 numSteps]);
    
    subplot(3,2,4);
    imagesc(1:numSteps, xbins(21:31), [postxM_r(21:30,:)+flipud(postxM_r(32:41,:)) ; ...
        postxM_r(31,:)]  );
    ylabel('Linearized position','FontSize',14);
    set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
    colormap(flipud(hot(256)));
    caxis([0 0.3]);
    xlim([0 numSteps]);
    
    subplot(3,2,6);
    imagesc(1:numSteps,xbins(ind_seg3),postxM_r(ind_seg3,:));
    ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
    set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
    colormap(flipud(hot(256)));
    caxis([0 0.3]);
    xlim([0 numSteps]);


elseif 1
    
    if eventoption == 2
        eventstring = 'ripplescons';
    elseif eventoption == 3;
        eventstring = 'wavegammacons';
    end
    
    xvecms = (1:numSteps)-ripms_start;
    ytetvec = 1:max(selected_tets);
    
    figure;
    
    % spike raster plot
    subplot(3,2,[1 3 5]);
    hold on
    for tet = 1:max(selected_tets)
        spikebins = find(spike_r(tet,:) > 0);
        numlines = length(spikebins);
        for ll = 1:numlines
            plot([spikebins(ll) spikebins(ll)] - ripms_start,[tet-0.5 tet+0.5],...
                'linestyle','-','Color',[0 0 0],'LineWidth',1)
        end
    end
    set(gca,'ytick',1:max(selected_tets))
    %set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
    xlabel('Time (ms)');ylabel('Tetrode');
    title([eventstring ' # : ',num2str(ripnum)],'fontweight','bold','fontsize',14);
    set(gca,'tickdir','out');

    box off
     % plot thick line indicating ripple
     hold on
    plot([ripms_start ripms_end] - ripms_start,[.5 .5],'Color','r','linewidth',4)
    if 0
        set(gca,'xtick',[-300:100:300]);
        set(gca,'xticklabel',{'-300','','','0','','','+300'});
        xlim([-300 300]);
    else
        set(gca,'xtick',[-500:100:500]);
        set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
        xlim([-500 500]);       
    end
    set(gca,'ydir','reverse')
    set(gca,'fontsize',12)    
    
    subplot(3,2,[2 4 6]);
    imagesc(xvecms,xbins,postxM_r);
    %title('postx, clusterless','FontSize',12);
    ylabel('Linearized position','FontSize',12);
    %set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
    %colormap(flipud(hot(256)));
    colormap(hot);
    caxis([0 0.2]);
    xlim([-300 300]);
    set(gca,'tickdir','out');
    if 0
        set(gca,'xtick',[-300:100:300]);
        set(gca,'xticklabel',{'-300','','','0','','','+300'});
        xlim([-300 300])
    else
        set(gca,'xtick',[-500:100:500]);
        set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
        xlim([-500 500])
    end
    set(gca,'fontsize',12)
    xlabel('Time (ms)')
    
    hold on
    % plot current position
    plot(rippos_timevec,rippos,'linewidth',2,'Color',[.6 .6 .6])
    
    % plot thick line indicating ripple
    plot([ripms_start ripms_end] - ripms_start,[-2.9 -2.9],'Color','r','linewidth',4)
    
    
end


pause
%close all
end
