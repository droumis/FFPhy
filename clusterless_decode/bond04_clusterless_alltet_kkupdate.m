
%%% Clusterless decoding script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%


datadir = '/opt/data50/kkay/DECODE/';

animalname = 'Bond';
day = 4;
ep = 2;
    animalinfo = animaldef(animalname);

SPIKETHRESH = 0;  % spikethresh (in uV, set to 0 if not want to use)
eventoption = 2;  %  ==0 for full epoch, == 2 for ripplescons, ==3 for WG
extra_time = 500;  % ms of extra time before and after to plot
ENCODETYPE = 2;    % == 1 for NAIVE, ==2 for STRICT

calculate = 0;
    outputfilename = 'Strict_0_clusterless';
    
plot_fullepoch = 1;
    if plot_fullepoch
        startsec = 525 ;        % time (in sec, within epoch) to begin plotting
        windowsec = 5 ;         % how many seconds to plot
        calc_events = 1;        % plots SWR and WG 
            timefilterscript
            eventspec_plot{1} = [ { rip2 } , {[1 .85 .85]}, {[1 .1 .1]} ] ;
            eventspec_plot{2} = [ { gf2_ca3dg }, {[.85 .85 1]}, {[.3 .3 1]} ];
        plot_powertrace = 1;    % plots ripple and WG power traces
    end

%%%%%% Specify event selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    timefilterscript

if eventoption == 0
    encodespec = { {norip, nogf_ca3dg, vel4}, [.7 .7 .7] };    
elseif eventoption == 2
    clear eventspec    
    eventspec = [ {rip6} , {[1 .6 .6]} , {[1 .1 .1]} ] ;   % {timefilter}, {patch color}, {line color}
    eventname = 'rip2';
    encodespec = {};
elseif eventoption == 3
    clear eventspec
    eventname = 'gammafcons';
    eventspec = [ {gf2_ca3dg} , {[.6 .6 1]}, {[.1 .1 1]} ] ;    
    encodespec = {};
end

if calculate

    tic
%%%%%% basic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',day);
linpos = loaddatastruct(animalinfo{2},animalinfo{3},'linpos',4);
postimevec=linpos{4}{2}.statematrix.time;            
poslin=vecLF(:,2);  
    poslin2 = poslin';

%%%%%% basic values:   xdel, xbins, dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%  Encode per tetrode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    inds2 = (filedata.params(:,2) > SPIKETHRESH) | (filedata.params(:,3) > SPIKETHRESH) | (filedata.params(:,4) > SPIKETHRESH) | (filedata.params(:,5) > SPIKETHRESH);
    % spikes that occur above 6 cm/s
    %spikeposinds = lookup(filedata.params(:,1),postimevec);
    %spikevels = pos{day}{ep}.data(spikeposinds,5);
    %inds3 = spikevels > 6;
    
    % determine final set of inds
    if ENCODETYPE == 1                                      % "NAIVE"
        inds = inds1 & inds2;   
        inds_encode = inds;
    elseif ENCODETYPE == 2                                    % STRICT
        timefilterscript;
        clear nonencodefilter
        nonencodefilter1 = {locomotion};
        nonencodefilter2 = {gf2_ca3dg};
        nonencodefilter3 = {rip2};
            % locomotion
        [EVENTPERIODS,TOTALTIME] = evaluatetimefilter(animalinfo{2}, animalinfo{3}, nonencodefilter1, [day ep]);
        nonencodeperiods1 = EVENTPERIODS{day}{ep};
        totaltime = TOTALTIME{day}{ep};
        %disp(sprintf('Filtering out of encoding: %d s locomotion filtered out',round(totaltime)));        
            % wave gamma
        [EVENTPERIODS,TOTALTIME] = evaluatetimefilter(animalinfo{2}, animalinfo{3}, nonencodefilter2, [day ep]);
        nonencodeperiods2 = EVENTPERIODS{day}{ep};
        totaltime = TOTALTIME{day}{ep};
        %disp(sprintf('Filtering out of encoding: %d s WG filtered out',round(totaltime)));                    
            % ripple
        [EVENTPERIODS,TOTALTIME] = evaluatetimefilter(animalinfo{2}, animalinfo{3}, nonencodefilter3, [day ep]);
        nonencodeperiods3 = EVENTPERIODS{day}{ep};
        totaltime = TOTALTIME{day}{ep};
        %disp(sprintf('Filtering out of encoding: %d s Ripple filtered out',round(totaltime)));                
        % spikes that occur during non-encoding periods
        inds4 = isExcluded(filedata.params(:,1)/10000,nonencodeperiods1);
        inds5 = isExcluded(filedata.params(:,1)/10000,nonencodeperiods2);
        inds6 = isExcluded(filedata.params(:,1)/10000,nonencodeperiods3);
        % intersection: valid spike inds
        inds = inds1 & inds2;
        inds_encode = inds1 & inds2 & ~inds4 & ~inds5 & ~inds6;
        disp(sprintf('%d of %d spikes excluded for being in non-encode period',sum(inds)-sum(inds_encode),sum(inds)))          
    end

    % mark vector (4 channel amplitudes)
    mark0{tet}=[filedata.params(inds,2) filedata.params(inds,3) filedata.params(inds,4) filedata.params(inds,5)];

    % spike timestamps
    spktimes{tet}=filedata.params(inds,1);   
    spktimes_encode{tet} = filedata.params(inds_encode,1)/10000; 
    
    [procInd0{tet},procInd1{tet}]=histc(spktimes_encode{tet},postimevec);   % for encoding
    [~            ,procInd1_ALL{tet}]=histc(spktimes{tet}/10000,postimevec);      % added this 6.27.15 kk to allow encoding spike set
    
    procInd{tet}=find(procInd0{tet});
    spikeT{tet}=postimevec(procInd{tet});
    spike{tet}=procInd0{tet}';
    [~,rawInd0{tet}]=histc(spktimes_encode{tet},spktimes_encode{tet});
    markAll{tet}(:,1)=procInd1{tet};
    markAll{tet}(:,2:5)=mark0{tet}(rawInd0{tet}(rawInd0{tet}~=0),:);
    mdel=20; ms = min(min(markAll{tet}(:,2:5))):mdel:max(max(markAll{tet}(:,2:5)));
    sxker=2*xdel; smker = mdel; T=size(postimevec,1);
    occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
    %occ: columns are identical; occupancy based on position; denominator
    Xnum{tet}=normpdf(xbins'*ones(1,length(spktimes_encode{tet})),ones(length(xbins),1)*poslin2(procInd1{tet}),sxker);
    %Xnum: Gaussian kernel estimators for position
    Lint{tet} = sum(Xnum{tet},2)./occ(:,1)./dt; %integral
    Lint{tet} = Lint{tet}./sum(Lint{tet});
    
end

% Collect all spikes, marks, and procInd, and lengths
spktimes_all = [];          % formerly time0 = [];   % --ALL-- spike times
spktimes_all_sorted = [];   % formerly time = [];    %   " " but sorted 
sinds = [];                 % formerly timeInd = [];  % all sorted spike inds
mark0_all = [];
procInd1_all = [];
len = {};
for tet = selected_tets
    
    numspktet = length(spktimes{tet});
    spktimes_all = [spktimes_all ;  spktimes{tet}    tet * ones(numspktet,1) ];   % * tag with the tet number
    mark0_all=[ mark0_all ; mark0{tet} ];
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
mdel=20; ms=100:mdel:max(mark0_all(:)); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_all=normpdf(xbins'*ones(1,length(spktimes_all)),ones(length(xbins),1)*poslin2(procInd1_all),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_all=sum(Xnum_all,2)./occ(:,1)./dt; %integral
Lint_all=Lint_all./sum(Lint_all);
%Lint: conditional intensity function for the unmarked case



%%%% Event identification (ripples or WG) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ti=round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector

if eventoption == 0
   
    
    plotwin_seg=[1 length(ti)]; %index for ripple segments

    extra_time = 0;
    
elseif eventoption == 2
    
    [EVENTPERIODS,TOTALTIME] = evaluatetimefilter(animalinfo{2}, animalinfo{3}, eventspec(1:(end-2)), [day ep]);
    riplist = EVENTPERIODS{day}{ep};
    totaltime = TOTALTIME{day}{ep};
        disp(sprintf('%d s of %s identified',round(totaltime),eventname));
    
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
    
    plotwin_seg=[find(startT)' find(endT)']; %index for ripple segments    
    
elseif eventoption == 3
    
  
    
end


%%
ripple_dur=plotwin_seg(:,2)-plotwin_seg(:,1);
%1203, 1237

%xi=round(time/10);  %spike times   % time is postimevec
% tot_ripple=zeros(length(ripple_dur),1);
% for i=1:length(ripple_dur)
%     numSteps=ripple_dur(i);
%     spike_tim=plotwin_seg(i,1):plotwin_seg(i,2); 
%     a1=0;
%     for t=1:numSteps
%         tt=spike_tim(t);
%         a=find(xi==ti(tt));
%         a1=a1+length(a);
%     end
%     tot_ripple(i)=a1;
% end


numplotwinseg = size(plotwin_seg,1);

for gg = 1:numplotwinseg
    
    % ms time bins to decode
    spike_tim = (plotwin_seg(gg,1)  - extra_time ) : (plotwin_seg(gg,2) + extra_time );      
    
    % ms time of ripple start
    ripms_start = lookup(plotwin_seg(gg,1),spike_tim);
    ripms_end = lookup(plotwin_seg(gg,2),spike_tim);
    
    % clock times for ripple
    rip_starttime = postimevec(1) + plotwin_seg(gg,1)/1000;
    rip_endtime = postimevec(1) + plotwin_seg(gg,2)/1000;
    
    % position of the animal during ripple (should barely change)
    if eventoption ~= 0
        PADDING = 20;  % how many position samples (29.97 Hz) before and after to collect pos for
    else
        PADDING = 0;
    end
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
            l0=normpdf(markAll{tet}(i,2)*ones(1,length(spktimes_encode{tet})),markAll{tet}(:,2)',smker)...   % is it supposed to be spktimes_encode length or spktimes ?
                .*normpdf(markAll{tet}(i,3)*ones(1,length(spktimes_encode{tet})),markAll{tet}(:,3)',smker)...
                .*normpdf(markAll{tet}(i,4)*ones(1,length(spktimes_encode{tet})),markAll{tet}(:,4)',smker)...
                .*normpdf(markAll{tet}(i,5)*ones(1,length(spktimes_encode{tet})),markAll{tet}(:,5)',smker);
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
    title(['ripple = ',num2str(gg)]);
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
    
    if eventoption == 0
        break   % don't print full epoch this way -- use plot_fullep code block
    elseif eventoption == 2
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
    title([eventstring ' # : ',num2str(gg)],'fontweight','bold','fontsize',14);
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
close all
end

    % if continuous epoch option, save
    if eventoption == 0
        cd(datadir)
        paraminfo = paramsstruct(encodespec);
        save(outputfilename,'spike_r','postxM_r','ti','xbins','vecLF','vecLF0','vecL0','paraminfo','-v7.3')
    end

    toc
end



if plot_fullepoch
    
    % basic information
    startsec;
    endsec = startsec + windowsec;
    tetinfo = loaddatastruct(animalinfo{2},animalinfo{3},'tetinfo');
    maxtet = max(selected_tets);
    
    % time variables
    epstart = ti(1);
    firstms = lookup(startsec*1000,ti-epstart);
    lastms = lookup(endsec*1000,ti-epstart);
    timevec_ms = ti(firstms:lastms);   % clock time vector
    
    H = figure('units','normalized','outerposition',[0 0 1 .4]);
    
    % Spike raster plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,2,[1 3 5]);
    hold on
    
    % plot SWR + WG as patches in background %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if calc_events
        eventperiods = {};
        goodinds = {};
        for evtype = 1:length(eventspec_plot)
            % grab the color to use
            eventclr = eventspec_plot{evtype}{2};
            % determine the event periods
            [EVENTPERIODS,TOTALTIME] = evaluatetimefilter(animalinfo{2}, animalinfo{3}, eventspec_plot{evtype}(1:(end-2)), [day ep]);
            eventperiods{evtype} = EVENTPERIODS{day}{ep};
            totaltime = TOTALTIME{day}{ep};
                disp(sprintf('%d s of event type #%d identified',round(totaltime),evtype));
            % identify which event periods overlap/are contained in this plotting window
            goodinds{evtype} = find( isExcluded(eventperiods{evtype}(:,1),timevec_ms/1000) | ...
                                     isExcluded(eventperiods{evtype}(:,2),timevec_ms/1000)   ); 
            numevents = length(goodinds{evtype});
            for ev = 1:numevents
                evind = goodinds{evtype}(ev);
                ev_winstart = eventperiods{evtype}(evind,1) - epstart/1000;
                ev_winend = eventperiods{evtype}(evind,2) - epstart/1000; 
                patch([ev_winstart ev_winstart ev_winend ev_winend],...
                      [0 maxtet+0.5 maxtet+0.5 0 ],eventclr,'edgecolor','none');
            end
        end
    else
        if ~exist('eventperiods','var')
            disp('eventperiods not detected')
        else
            disp('using previously calculated eventperiods')
        end
    end
    
    % plot rip and wg powertrace in background %%%%%%%%%%%%%%%%%%%%%%%%
    if plot_powertrace
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', day, ep);
            riptrace = zscorer(out{day}{ep}{1}.powertrace);
            	tvec_rip = out{day}{ep}{1}.eegtimesvec_ref;
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', day, ep);
            wgtrace = zscorer(out{day}{ep}{2}.powertrace);
            	tvec_wg = out{day}{ep}{2}.eegtimesvec_ref;
        if ~all(tvec_wg == tvec_rip)
            error('these need to be the same')
        end
        indA = lookup(startsec + epstart/1000,tvec_rip);
        indB = lookup(endsec + epstart/1000,tvec_rip);
        % plot rip
        plot(tvec_rip(indA:indB)-epstart/1000,...
             -riptrace(indA:indB)  + maxtet + 4,...
                'linewidth',2,'Color',[1 .3 .3]);
        % plot WG
         plot(tvec_wg(indA:indB)-epstart/1000,...
              -wgtrace(indA:indB) + maxtet + 4,...
                'linewidth',2,'Color',[.3 .3 1]);       
    end
    
    % Plot spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for tet = 1:maxtet
        % identify spikes in window
        %spikebins = find(spike_r(tet,firstms:lastms) > 0);
        spikebins = find(spike_r(tet,firstms:lastms) > 0) / 1000 + startsec;

        numlines = length(spikebins);
        if isfield(tetinfo{day}{ep}{tet},'area')
            if strcmp(tetinfo{day}{ep}{tet}.area,'CA1')
                spkclr = [0 0 0];
            elseif strcmp(tetinfo{day}{ep}{tet}.area,'CA3')
                spkclr = [1 .3 .3];
            end
        else
            spkclr = [.7 .7 .7];
        end
        for ll = 1:numlines
            plot([spikebins(ll) spikebins(ll)],[tet-0.5 tet+0.5],...
                'linestyle','-','Color',spkclr,'LineWidth',1)
        end
    end
    
    % formatting the plot %%%%%%%%%%%
    ylim([0 maxtet+0.5 + 5]);
    xlim([startsec  endsec]);    
    set(gca,'ytick',1:maxtet)
    set(gca,'xtick',[startsec:1:endsec]);        
    title('Continuous decoding','fontweight','bold','fontsize',14);
    set(gca,'tickdir','out');
    xlabel('Time (ms)');    ylabel('Tetrode');       
    box off
    set(gca,'ydir','reverse')
    set(gca,'fontsize',12)    
    
    

    % Decoding plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,2,[2 4 6]);
    
    % plot posteriors
    imagesc(timevec_ms/1000 - epstart/1000, xbins, postxM_r(:,firstms:lastms));  hold on;
    
    % formatting
    ylabel('Linearized position','FontSize',12);
    xlabel('Time (ms)')
    colormap(hot);
    caxis([0 0.2]);
    xlim([startsec  endsec]);
    set(gca,'xtick',[startsec:1:endsec]);    
    set(gca,'tickdir','out');
    set(gca,'fontsize',12)
    
    % plot current position
    posA = lookup(startsec + epstart/1000,vecLF(:,1));  % first column is timevec (postimevec)
    posB = lookup(endsec + epstart/1000,vecLF(:,1));
    
    plot( vecLF(posA:posB,1) - epstart/1000, vecLF(posA:posB,2),...
          'linewidth',2,'Color',[.6 .6 .6])
    
    % Plot events as thick lines (uses data from calc_events above)
    if 1
        eventperiods;
        goodinds;
        for evtype = 1:length(eventspec_plot)
            % grab the color to use
            eventclr = eventspec_plot{evtype}{3};
            numevents = length(goodinds{evtype});
            for ev = 1:numevents
                evind = goodinds{evtype}(ev);
                ev_winstart = eventperiods{evtype}(evind,1) - epstart/1000;
                ev_winend = eventperiods{evtype}(evind,2) - epstart/1000;
                % 
                if 1
                    plot([ev_winstart ev_winend],[-2.9 -2.9],'Color',eventclr,'linewidth',4)
                    hold on
                else
                    patch([ev_winstart ev_winstart ev_winend ev_winend],...
                          [-3 3 3 -3 ],eventclr,'edgecolor','none');
                end
            end
        end
    end
    
end






















