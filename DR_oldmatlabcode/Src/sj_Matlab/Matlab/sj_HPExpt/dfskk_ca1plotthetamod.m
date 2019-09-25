
% From dfskk_ca1plotthetamod
% Plot theta modulation of a population of cells - either PFC or CA1

%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Ten'};
%animals = {'Dudley'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Frank'};
animals = {'Egypt'};
%animals = {'Bond','Corriander','Frank','Egypt','Chapati'};

%-----------------------------------------------------

runscript = 0;
%statictet = 11;   % theta eeg reference for all specified days
     % Chapati: best staticref is tetrode 1. matches Mizuseki for CA1
     % Egypt: best staticref is tetrode 11. matches Mizuseki for CA1
     % Frank: tetrode 30
     % Bond: tetrode 30
     % Corriander: reference is tetrode 24
     % Barack: reference is 4 for days 1-17
     % Calvin: reference is 16
     % Dwight; reference is 7

 
if runscript == 1
    
    clear eegfilter
    clear super;
%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = '(isequal($type, ''run'') || isequal($type,''sleep''))';
%epochfilter{1} = '(($exposure > 0 ))';
%epochfilter{1} = '(isequal($type, ''run''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 7))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 7))';

ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 7))';
%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

%timefilter_ca1 = { {'get2dstate', '((abs($velocity) >= 8))'} };
%timefilter_ca3 = { {'get2dstate', '((abs($velocity) >= 8))'} };
timefilter_ca1 = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                  {'kk_get2dstate', '((abs($velocity) > 5))'}, ... 
                  {'gethighthetatimes2','($ntheta >= 2)',[],'tetfilter','(isequal($area, ''CA1'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};
            % {'gethighthetatimes2','($ntheta == 1)',30,'powerratio1',4,'powerratio2',.3,'mindur',.5}};
            % {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}, ...
timefilter_ca3 = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                  {'kk_get2dstate', '((abs($velocity) > 5))'}, ... 
                  {'gethighthetatimes2','($ntheta >= 2)',[],'tetfilter','(isequal($area, ''CA3'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};
timefilter_ca2 = {{'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                  {'kk_get2dstate', '((abs($velocity) > 5))'}, ... 
                  {'gethighthetatimes2','($ntheta >= 2)',[],'tetfilter','(isequal($area, ''CA1'') && ($numcells >= 1))','powerratio1',4,'powerratio2',.3,'mindur',.5}};

%eegfilter = {'geteegtet', 'theta', 'statictet',statictet};
%ca1tetfilter = '(isequal($area, ''CA1''))';
%eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter',ca1tetfilter};
                   % picks tetrode with largest theta variance 
%eegfilter = {'geteegtet', 'theta', 'sametet', 1};
%eegfilter = {'geteegtet', 'theta', 'statictet',statictet};
eegfilter = {'geteegtet', 'thetagnd', 'sametet',1}

iterator = 'singlecelleeganal';

ca1f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter_ca1, 'cells', ca1cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter_ca3, 'cells', ca3cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
ca2f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter_ca2, 'cells', ca2cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);

ca1f = setfilterfunction(ca1f, 'dfakk_thetagndphase', {'spikes', 'thetagnd'},'nbins',24);
ca3f = setfilterfunction(ca3f, 'dfakk_thetagndphase', {'spikes', 'thetagnd'},'nbins',24);
ca2f = setfilterfunction(ca2f, 'dfakk_thetagndphase', {'spikes', 'thetagnd'},'nbins',24);
%dummy = setfilterfunction(ca1f, 'plotthetamod', {'spikes', 'theta'},'color','k','nbins',36);
%ca2f = setfilterfunction(ca2f, 'dfakk_thetagndphase', {'spikes', 'theta'},'color','g','nbins',36);

ca2f = runfilter(ca2f);
ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

end

%% Designate hippocampal area to analyze.

for reg=[2]
    if reg==1
        region='CA1';
        clr=[0.3 0.3 0.3];
        caf=ca1f;
    elseif reg==2
        region='CA2';
        clr='g';
        caf=ca2f;
    else
        region='CA3';
        clr='r';
        caf=ca3f;
    end

%% Generate reference string for labelling figures.

if isequal(eegfilter,{'geteegtet', 'theta', 'sametet', 1})      % local reference
    referencestring = 'sametet';
elseif isequal(eegfilter,{'geteegtet', 'thetagnd', 'sametet', 1})
    referencestring = sprintf('sametet, ref to gnd');     
elseif isequal(eegfilter,{'geteegtet', 'theta', 'statictet', statictet})
    referencestring = sprintf('static ref, tetrode %d',statictet);       % static reference
elseif isequal(eegfilter,{'geteegtet', 'thetagnd', 'statictet', statictet})    
    referencestring = sprintf('static ref, tetrode %d',statictet);       % static reference
end
    
%% Consolidate single cells' across epochs in a day (.celloutput field)

dummyindex=[];
for i=1:length(caf.output{1})
    dummyindex = [dummyindex;caf.output{1}(i).index];           % collect all day-epoch indices
end
caf.celloutput=struct;

for i=1:size(caf.output{1},2)
    daytetcell=caf.output{1}(i).index([1 3 4]);
    ind=[];
    while rowfind(daytetcell,dummyindex(:,[1 3 4]))~=0          % collect all rows (epochs)
        ind = [ind rowfind(daytetcell,dummyindex(:,[1 3 4]))];
        dummyindex(rowfind(daytetcell,dummyindex(:,[1 3 4])),:)=[0 0 0 0];
    end
     
    dummy=[];
    for r=ind
        dummy=[dummy ; caf.output{1}(r).sph];
    end

    if ~isempty(dummy)
        if ~isfield(caf.celloutput,'index')      % if initial entry
            caf.celloutput(1).index=daytetcell;
            caf.celloutput(1).spikephases=dummy;
        else            
            A.index=daytetcell;
            A.spikephases=dummy;
            caf.celloutput=horzcat(caf.celloutput,A);
        end
    end
    
end

%% Compute mean phase, modulation depth, phase histograms, Rayleigh test, and plot overall modulation w/ SE.

% set # bins
nbins=24;       
    bins = -pi:(2*pi/nbins):pi;
% compute mean phase and phase histograms
clear i;
for k=1:length(caf.celloutput)
    A=mean(exp(i*caf.celloutput(k).spikephases));
    meancos=real(A);
    meansin=imag(A);
    caf.celloutput(k).meanphase=atan2(meansin,meancos);
    caf.celloutput(k).moddepth=abs(A);
    caf.celloutput(k).bins=bins;
    caf.celloutput(k).phasehist=histc(caf.celloutput(k).spikephases,bins);
    out=rayleigh_test(caf.celloutput(k).spikephases);
    caf.celloutput(k).rayleighp=out.p;
end

    %% add to "super" variable
    
super{reg}=caf;



%% Plots.


if 1
%plot phase histogram of mean phases
figure
dummy=[]; sigcount=0;
    %collect all units' mean phases
for k=1:length(caf.celloutput)
    if (caf.celloutput(k).rayleighp < 0.05)
        dummy=[dummy ; caf.celloutput(k).meanphase];
        sigcount=sigcount+1;
    end
end
sigcount
N = histc(dummy,bins);
    N=N(1:(end-1));
    bins_plot = bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
h = bar([bins_plot bins_plot+2*pi],[N ; N],'histc');
    set(h(1),'facecolor',clr)
    set(h(1),'edgecolor',clr)
title(sprintf('%s %s mean phases of all units // %s',animals{1},region,referencestring),'FontSize',16,'FontWeight','bold')
axis tight
hold on
plot([pi,pi],[0 999],'k--','LineWidth',1.5)
plot([-pi,-pi],[0 999],'k--','LineWidth',1.5)
plot([3*pi,3*pi],[0 999],'k--','LineWidth',1.5)
plot([0,0],[0 999],'k:','LineWidth',1.5)
plot([2*pi,2*pi],[0 999],'k:','LineWidth',1.5)
end


if 0
%plot distribution of modulation depths
figure
moddepths=[];
for i=1:length(caf.celloutput)
   moddepths=[moddepths ; caf.celloutput(i).moddepth];
end
hist(moddepths,10,clr)
title(sprintf('%s distribution of modulation depths',region))
end




if 0
% plot phase histogram of aggregate spikes, sig units

norm = 1;

figure
dummy=[]; sigcount=0;
for k=1:length(caf.celloutput)
    if 1 % (caf.celloutput(k).rayleighp < 0.05)       % check if significant
        dummy=[dummy ; caf.celloutput(k).spikephases];
        sigcount=sigcount+1;
    end
end
sigcount
phasehist=histc(dummy,bins);
    phasehist = phasehist(1:(end-1));
        phasehist_norm = phasehist/(sum(phasehist));
    bins_plot = bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
    if norm == 1
        h = bar([bins_plot bins_plot+2*pi],[phasehist_norm; phasehist_norm],'histc');
        set(h(1),'facecolor',clr)
        set(h(1),'edgecolor',clr)
        axis tight
        ylim([0 .1])        
    else
        h = bar([bins_plot bins_plot+2*pi],[phasehist; phasehist],'histc');
        set(h(1),'facecolor',clr)
        set(h(1),'edgecolor',clr)
        axis tight
        ylim([0 max(phasehist)+1500])
    end
titlestring=sprintf('%s %s phase hist of all spikes, sig units // %s',animals{1},region,referencestring);
title(titlestring,'FontSize',14,'FontWeight','bold')
hold on
plot([pi,pi],[0 99999],'k--','LineWidth',1.5)
plot([-pi,-pi],[0 99999],'k--','LineWidth',1.5)
plot([3*pi,3*pi],[0 99999],'k--','LineWidth',1.5)
plot([0,0],[0 99999],'k:','LineWidth',1.5)
plot([2*pi,2*pi],[0 99999],'k:','LineWidth',1.5)

end


if 1
% plot normalized histogram sig units, SEM

figure
dummy=[]; sigcount=0;
for k=1:length(caf.celloutput)
    if (caf.celloutput(k).rayleighp < 0.05)       % check if significant
        dummy=[dummy ; caf.celloutput(k).phasehist'/sum(caf.celloutput(k).phasehist)];   %normalize
        sigcount=sigcount+1;
    end
end
sigcount
bins = caf.celloutput(k).bins;
phasehist_mean=mean(dummy,1);
phasehist_sem=std(dummy,1)/sqrt(size(dummy,1));
    phasehist_mean = phasehist_mean(1:(end-1));
    bins_plot = bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
       
    h = bar([bins_plot bins_plot+2*pi],[phasehist_mean phasehist_mean],'histc');
    set(h(1),'facecolor',clr)
    set(h(1),'edgecolor',clr)
    axis tight
    ylim([0 .1])
    hold on
    
    % plot sem bars
    for jj=1:length(bins_plot)
        plot([bins_plot(jj),bins_plot(jj)],[phasehist_mean(jj)-phasehist_sem(jj) phasehist_mean(jj)+phasehist_sem(jj)],'k','LineWidth',1.5)
        plot([bins_plot(jj)+2*pi,bins_plot(jj)+2*pi],[phasehist_mean(jj)-phasehist_sem(jj) phasehist_mean(jj)+phasehist_sem(jj)],'k','LineWidth',1.5)
    end

titlestring=sprintf('%s %s phase hist of all spikes, sig units // %s',animals{1},region,referencestring);
title(titlestring,'FontSize',14,'FontWeight','bold')
hold on
plot([pi,pi],[0 99999],'k--','LineWidth',1.5)
plot([-pi,-pi],[0 99999],'k--','LineWidth',1.5)
plot([3*pi,3*pi],[0 99999],'k--','LineWidth',1.5)
plot([0,0],[0 99999],'k:','LineWidth',1.5)
plot([2*pi,2*pi],[0 99999],'k:','LineWidth',1.5)

end



if 1
% plot individual phase histogram of all units

norm = 1;

figure
titlestring=sprintf('%s %s phase hist of individual units // %s',animals{1},region,referencestring);
title(titlestring,'FontSize',14,'FontWeight','bold')
counter=1;
for k=1:length(caf.celloutput)
    if counter==81
        counter=1;
        figure
        titlestring=sprintf('%s %s phase hist of individual units // %s',animals{1},region,referencestring);
        title(titlestring,'FontSize',14,'FontWeight','bold')
    end
    subplot(8,10,counter)
    bins_plot = caf.celloutput(k).bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
    phasehist = caf.celloutput(k).phasehist(1:(end-1));
        phasehist_norm = phasehist/sum(phasehist);
    if norm == 1
        if size(phasehist_norm,1) < size(phasehist_norm,2)
            phasehist_norm = phasehist_norm';
        end
            %plot
            h = bar([bins_plot bins_plot+2*pi],[phasehist_norm ; phasehist_norm],'histc');
            title(num2str(caf.celloutput(k).index))
            axis tight
            ylim([0 .2])
    else
        if size(phasehist,1) < size(phasehist,2)
            phasehist = phasehist';
        end
            %plot
            h = bar([bins_plot bins_plot+2*pi],[phasehist ; phasehist],'histc');
            title(num2str(caf.celloutput(k).index),'FontSize',12,'FontWeight','bold')
            axis tight
            ylim([0 250])
    end
    
    set(h(1),'facecolor',clr)
    set(h(1),'edgecolor',clr)

    % plot guide lines
    hold on
    plot([pi,pi],[0 9999],'k--','LineWidth',1.5)
    plot([-pi,-pi],[0 9999],'k--','LineWidth',1.5)
    plot([3*pi,3*pi],[0 9999],'k--','LineWidth',1.5)
    plot([0,0],[0 9999],'k:','LineWidth',1.5)
    plot([2*pi,2*pi],[0 9999],'k:','LineWidth',1.5)
    
    counter=counter+1;
end
end



end


% % bar
% count = histc(allspikephases, bins);
% out = bar(bins, count, 'hist');
% set(out,'facecolor','k')
% title('aggregate theta modulation');
% 
% % lineplot
% dischargeprob=count./sum(count);
% plot(bins(1:(end-1)),dischargeprob(1:(end-1)),'k','LineWidth',2);
% 
% [m ph] = modulation(allspikephases);






