
% plots aggregate spiking relative to gamma cycle peaks nearest gamma power
% peaks of extracted gamma "events" (2 SD)

runscript = 0

% pre parameters
animal_list = {'Chapati'}
minthresh = 0;     % how big the gammah events are, in SD (usually 2)
nbins = 24;
ngammah = 0;       % usually 1, if not more
eventreg = 'CA1';
phasereg = 'CA1';
    % phase reference approaches
ground = 1;
onephasetet = 0;   % 1. if onephasetet and localtet = 0, then
localtet = 1;      % dfakk_getgammatrig spiking detects CA1 gamma events
                   % from ANY CA1 tetrode, and chooses the CA1 tetrode with
                   % the most cells + high power participation in the event
                   % 2. if onephasetet = 1, then tells kk_getgammahtimes to pick one tetrode (most # cells) per day for events AND phase. For this, must set ngammah to 1.
                   % 3. if localtet = 1, this tells dfakk_getgammatrigspiking to
                                %    i.   only count spikes from locally detected gammah
                                %    ii.  only use phase from local tetrode gamma
                                %    iii. for CA3, only consult the CA1
                                %    tetrode with which it has the highest coherence for
                                %    that epoch

epoch = 'all';
timefilt_flag = '11';  % first digit is ripple exclusion (1) or not (0)
                        % second digit is > velocity threshold   -- ex. 00 is most inclusive
% analyses
individual_flag = 0;
aggregate_flag = 1;

% post parameters
regions = [1 2 3];
normalize = 1;
sig_only = 1;
minimum_spikes = 10;
ca2_unitlist = [4 3 2 ; 4 17 2 ; 4 17 5 ; 5 17 2; 6 17 4 ; 7 17 5 ; 8 3 1; 9 17 5 ; 10 3 1; 10 3 3; 10 3 4; 10 3 7; 10 3 10];
ca2_classification_flag = 1;
            
if runscript

% Animal Selection
animals = animal_list;

% Epoch Filter
if strcmp(epoch,'run')
    epochfilter = '(isequal($type, ''run''))';
elseif strcmp(epoch,'sleep')
    epochfilter = '(isequal($type, ''sleep''))';
elseif strcmp(epoch,'all')
    epochfilter = '(isequal($type, ''run'') || isequal($type, ''sleep''))';
else
    error('specify epoch type!')
end

% Time Filter
if strcmp(timefilt_flag,'00')         % all times
    timefilter = { {'kk_get2dstate', '($velocity >= 0)'} };
elseif strcmp(timefilt_flag,'10')     % all times but ripples excluded
    timefilter = { {'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'} };
elseif strcmp(timefilt_flag,'11')     % exploratory periods
    timefilter = { {'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                   {'kk_get2dstate', '($velocity >= 4)'} };
elseif strcmp(timefilt_flag,'rip')    % ripple periods
    timefilter = { {'kk_getriptimes', '($nripples >= 2)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                    {'kk_get2dstate', '($velocity < 4)'}  };                
end
    
% 
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';    %% ($meanrate < 3))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';   
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 3))';
%ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && ($meanrate < 3))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 3))';

% Iterator
iterator = 'kk_multicellanal';

% Filter Creation
ca1f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca2f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
if ground
ca1f = setfilterfunction(ca1f, 'dfakk_gethighgammagndtrigspiking2', {'highgammagnd','gammah','spikes','tetinfo'},'gammatetfilter',['(isequal($area, ''' phasereg ''') && ($numcells >= 2))'],'local',1);
ca2f = setfilterfunction(ca2f, 'dfakk_gethighgammagndtrigspiking2', {'highgammagnd','gammah','spikes','tetinfo'},'gammatetfilter',['(isequal($area, ''' phasereg ''') && ($numcells >= 2))'],'local',1);
ca3f = setfilterfunction(ca3f, 'dfakk_gethighgammagndtrigspiking2', {'highgammagnd','gammah','spikes','tetinfo'},'gammatetfilter',['(isequal($area, ''' phasereg ''') && ($numcells >= 2))'],'local',1);
else
ca1f = setfilterfunction(ca1f, 'dfakk_gethighgammatrigspiking2', {'highgamma','gammah','spikes','tetinfo'},'gammatetfilter',['(isequal($area, ''' phasereg ''') && ($numcells >= 2))'],'local',1);
ca2f = setfilterfunction(ca2f, 'dfakk_gethighgammatrigspiking2', {'highgamma','gammah','spikes','tetinfo'},'gammatetfilter',['(isequal($area, ''' phasereg ''') && ($numcells >= 2))'],'local',1);
ca3f = setfilterfunction(ca3f, 'dfakk_gethighgammatrigspiking2', {'highgamma','gammah','spikes','tetinfo'},'gammatetfilter',['(isequal($area, ''' phasereg ''') && ($numcells >= 2))'],'local',1);
end

% Run Analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end



% Parameter string (for figures).
params_string = [epoch ' epochs / ' 'timefilter: ' timefilt_flag ' / minthresh: ' num2str(minthresh) ...
                 ' / ngammah: ' num2str(ngammah) ' / phaseregion: ' phasereg ' / onephasetet: ' num2str(onephasetet)]
               

             

% Consolidate unit data over epochs

for reg=regions       % iterate over regions
    if reg==1
        region='CA1';
        clr=[0 0 0];
        f=ca1f;
    elseif reg==2
        region='CA2';
        clr=[0 1 0];
        f=ca2f;
    else
        region='CA3';
        clr=[1 0 0];
        f=ca3f;
    end

% Consolidate cells across epochs, within a day

%collect all cell indices, ignoring specific epochs, into daytetcell
indices = [];
for i=1:length(f.output{1})  % iterate over epochs
    indices = [indices ; f.output{1}(i).indices_coherence];
end
daytetcell = unique(indices(:,[1 3 4]),'rows');

f.celloutput = struct;

for ind=1:size(daytetcell,1)
    f.celloutput(ind).daytetcell = daytetcell(ind,:);
    f.celloutput(ind).spikephases = [];
    f.celloutput(ind).numspikes_total = 0;
    f.celloutput(ind).numspikes_gammah = 0;
    %consolidate data across epochs
    for ep=1:length(f.output{1})
        if ~isempty(f.output{1}(ep).indices_coherence)   % epoch may not have any data (fully excluded)
            jj=rowfind(daytetcell(ind,:),f.output{1}(ep).indices_coherence(:,[1 3 4]));
            if jj
                f.celloutput(ind).spikephases = [f.celloutput(ind).spikephases ; f.output{1}(ep).spikephases{jj}];
                f.celloutput(ind).numspikes_total = f.celloutput(ind).numspikes_total + f.output{1}(ep).numspikes_total(jj);
                f.celloutput(ind).numspikes_gammah = f.celloutput(ind).numspikes_gammah + f.output{1}(ep).numspikes_gammah(jj);
            end
        end
    end
end
    
    super{reg} = f;
   
end



%% Compute mean phase, modulation depth, phase histograms, Rayleigh test, and plot overall modulation w/ SE.

% set # bins

bins = -pi:(2*pi/nbins):pi;

% compute mean phase and phase histograms
clear i;
for reg=regions
    for k=1:length(super{reg}.celloutput)
        A=mean(exp(i*double(super{reg}.celloutput(k).spikephases)));
        meancos=real(A);
        meansin=imag(A);
        super{reg}.celloutput(k).meanphase=atan2(meansin,meancos);
        super{reg}.celloutput(k).moddepth=abs(A);
        super{reg}.celloutput(k).bins=bins;
        super{reg}.celloutput(k).phasehist=histc(double(super{reg}.celloutput(k).spikephases) / 10000,bins);
        out=rayleigh_test(double(super{reg}.celloutput(k).spikephases / 10000));
        super{reg}.celloutput(k).rayleighp=out.p;
    end
end




%% Plots.

if individual_flag
    
% plot individual phase histogram of all units

for reg = regions
    
    if reg==1
        region='CA1';
        clr=[0 0 0];
    elseif reg==2
        region='CA2';
        clr=[0 1 0];
    else
        region='CA3';
        clr=[1 0 0];
    end
    

figure
titlestring=sprintf('%s %s phase hist of individual units',animals{1},region);
title({titlestring params_string},'FontSize',14,'FontWeight','bold')
counter=1;
for k=1:length(super{reg}.celloutput)
    
    % don't plot if too few spikes
    if super{reg}.celloutput(k).numspikes_gammah < minimum_spikes
        continue
    end
    
    if counter==81
        counter=1;
        figure
        titlestring=sprintf('%s %s phase hist of individual units',animals{1},region);
        title({titlestring params_string},'FontSize',14,'FontWeight','bold')
    end
    subplot(8,10,counter)
    bins_plot = super{reg}.celloutput(k).bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
    phasehist = super{reg}.celloutput(k).phasehist(1:(end-1));
        phasehist_norm = phasehist/sum(phasehist);
        
    if (sig_only == 1) && (super{reg}.celloutput(k).rayleighp > 0.05)
        continue
    end
        
    if normalize == 1
        if size(phasehist_norm,1) < size(phasehist_norm,2)
            phasehist_norm = phasehist_norm';
        end
            %plot
            if super{reg}.celloutput(k).rayleighp <= 0.05
                fontcolor = [1 .1 .1];
                disp('significant, color')
            else
                fontcolor = [0 0 0]; 
            end
            h = bar([bins_plot bins_plot+2*pi],[phasehist_norm ; phasehist_norm],'histc');
            title(num2str(super{reg}.celloutput(k).daytetcell),'FontSize',12,'FontWeight','bold','Color',fontcolor)
            axis tight
            %ylim([0 .2])
    else
        if size(phasehist,1) < size(phasehist,2)
            phasehist = phasehist';
        end
            %plot
            h = bar([bins_plot bins_plot+2*pi],[phasehist ; phasehist],'histc');
            if super{reg}.celloutput(k).rayleighp <= 0.05
                fontcolor = [1 .1 .1];
                disp('significant, color')
            else
                fontcolor = [0 0 0]; 
            end
            title(num2str(super{reg}.celloutput(k).daytetcell),'FontSize',12,'FontWeight','bold','Color',fontcolor)
            axis tight
            %ylim([0 250])
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


if aggregate_flag
% plot phase histogram of aggregate spikes, sig units
for reg = regions
   

    if reg==1
        region='CA1';
        clr=[0 0 0];
    elseif reg==2
        region='CA2';
        clr=[0 1 0];
    else
        region='CA3';
        clr=[1 0 0];
    end


figure

dummy=[]; dummy2=[]; sigcount=0; cellcount = 0;

% collect spikes
for k=1:length(super{reg}.celloutput)
    if sig_only
        % check if significant and at least minimum # of spikes
        if (super{reg}.celloutput(k).rayleighp <= 0.05) && (super{reg}.celloutput(k).numspikes_gammah >= minimum_spikes)
            dummy=[dummy ; double(super{reg}.celloutput(k).spikephases) / 10000];
            sigcount=sigcount+1;
        end
    else
        dummy=[dummy ; double(super{reg}.celloutput(k).spikephases) / 10000];
        cellcount=cellcount+1;       
    end
end

cellcount


    % if normalize then can plot SEM
    if normalize == 1
        for k=1:length(super{reg}.celloutput)
            if ca2_classification_flag
            if strcmp(region,'CA2') && ~isempty(ca2_unitlist)
                if rowfind(super{reg}.celloutput(k).daytetcell,ca2_unitlist)
                    continue
                end
            end
            end
            if sig_only
                if (super{reg}.celloutput(k).rayleighp < 0.05)       % check if significant
                    dummy2=[dummy2 ; super{reg}.celloutput(k).phasehist'/sum(super{reg}.celloutput(k).phasehist)];   %normalize
                    sigcount=sigcount+1;
                end
            else
                dummy2=[dummy2 ; super{reg}.celloutput(k).phasehist'/sum(super{reg}.celloutput(k).phasehist)];   %normalize
            end

        end
        
        % truncate last bins
        dummy2 = dummy2(:,1:(end-1));
        % SEM
        phasehist_mean=mean(dummy2,1);
        phasehist_sem=std(dummy2,1)/sqrt(size(dummy2,1));
        bins_plot = bins(1:(end-1));
        bins_plot = bins_plot + (bins(2)-bins(1))/2;
        
        h = bar([bins_plot bins_plot+2*pi],[phasehist_mean phasehist_mean],'histc');
        hold on
        set(h(1),'facecolor',clr)
        set(h(1),'edgecolor',clr)
        axis tight
        ylim([0 .2])
        % plot sem bars
        for jj=1:length(bins_plot)
            plot([bins_plot(jj),bins_plot(jj)],[phasehist_mean(jj)-phasehist_sem(jj) phasehist_mean(jj)+phasehist_sem(jj)],'k','LineWidth',1.5)
            plot([bins_plot(jj)+2*pi,bins_plot(jj)+2*pi],[phasehist_mean(jj)-phasehist_sem(jj) phasehist_mean(jj)+phasehist_sem(jj)],'k','LineWidth',1.5)
        end
    % otherwise pile all spikes into the histogram
    else
        phasehist=histc(double(dummy(:)),bins);
        phasehist = phasehist(1:(end-1));
        bins_plot = bins(1:(end-1));
        bins_plot = bins_plot + (bins(2)-bins(1))/2;
        % plot
        h = bar([bins_plot bins_plot+2*pi],[phasehist; phasehist],'histc');
        set(h(1),'facecolor',clr)
        set(h(1),'edgecolor',clr)
        axis tight
        ylim([0 max(phasehist)+100])
    end
titlestring=sprintf('%s %s phase hist of all spikes, sig units',animals{1},region);
title({titlestring params_string},'FontSize',14,'FontWeight','bold')
hold on
plot([pi,pi],[0 99999],'k--','LineWidth',1.5)
plot([-pi,-pi],[0 99999],'k--','LineWidth',1.5)
plot([3*pi,3*pi],[0 99999],'k--','LineWidth',1.5)
plot([0,0],[0 99999],'k:','LineWidth',1.5)
plot([2*pi,2*pi],[0 99999],'k:','LineWidth',1.5)
end
end






if 0
    
for reg = regions;    
    
%plot phase histogram of mean phases
figure

dummy=[]; sigcount=0;
    %collect all units' mean phases
for k=1:length(super{reg}.celloutput)
    if (super{reg}.celloutput(k).rayleighp < 0.05)
        dummy=[dummy ; super{reg}.celloutput(k).meanphase];
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
title(sprintf('%s %s mean phases of all units',animals{1},region),'FontSize',16,'FontWeight','bold')
axis tight
hold on
plot([pi,pi],[0 999],'k--','LineWidth',1.5)
plot([-pi,-pi],[0 999],'k--','LineWidth',1.5)
plot([3*pi,3*pi],[0 999],'k--','LineWidth',1.5)
plot([0,0],[0 999],'k:','LineWidth',1.5)
plot([2*pi,2*pi],[0 999],'k:','LineWidth',1.5)
end
end


if 0
%plot distribution of modulation depths
for reg = regions
figure
moddepths=[];
for i=1:length(super{reg}.celloutput)
   moddepths=[moddepths ; super{reg}.celloutput(i).moddepth];
end
hist(moddepths,10,clr)
title(sprintf('%s distribution of modulation depths',region))
end
end





if 0
% plot normalized histogram sig units, SEM

figure
dummy=[]; sigcount=0;
for k=1:length(super{reg}.celloutput)
    if (super{reg}.celloutput(k).rayleighp < 0.05)       % check if significant
        dummy=[dummy ; super{reg}.celloutput(k).phasehist'/sum(super{reg}.celloutput(k).phasehist)];   %normalize
        sigcount=sigcount+1;
    end
end
sigcount
bins = super{reg}.celloutput(k).bins;
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

titlestring=sprintf('%s %s phase hist of all spikes, sig units',animals{1},region);
title({titlestring params_string},'FontSize',14,'FontWeight','bold')
hold on
plot([pi,pi],[0 99999],'k--','LineWidth',1.5)
plot([-pi,-pi],[0 99999],'k--','LineWidth',1.5)
plot([3*pi,3*pi],[0 99999],'k--','LineWidth',1.5)
plot([0,0],[0 99999],'k:','LineWidth',1.5)
plot([2*pi,2*pi],[0 99999],'k:','LineWidth',1.5)

end



















