
datadir = '/opt/data13/kkay/Supertheta_data/';


setparams = 1;
runscript = 0;
postprocessing = 0;  % consolidate over epochs to get .celloutput + perform calculations

% plots
plot_results = 1;
    if plot_results
        plot_singles = 0;
        plot_population = 1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if setparams
    nbins = 24;
        bins = -pi:(2*pi/nbins):pi;
    % define states
    timefilterscript; clear statespec;       
    statespec{1} = {{ fullepoch, norip }, 'runlinear_all'};            % all times          
    statespec{2} = {{ locomotion }, 'runlinear_all'} ;      % moving
    statespec{3} = {{ highspeed }, 'runlinear_all'} ;          
    statespec{4} = {{ immobile2, norip  }, 'runlinear_all'} ; 
    statespec{5} = {{ ripvel, norip  }, 'runlinear_all'} ; 
    statespec{6} = {{ locomotion  }, 'sleep_inclusive'} ; 
    statespec{7} = {{ highspeed  },'sleep_inclusive'} ; 
    statespec{8} = {{ sleepc  }, 'sleep_inclusive'} ; 
    statespec{9} = {{ sleepc,norip  }, 'sleep_inclusive'} ; 
    
    numstates = length(statespec);
    % eeg selection  (choose one or the other)
    sametet = 0;
    statictet = 1;   % theta eeg reference for all specified days - the selection for each animal is in the runscript section
    statictet_region = 'CA1';
    if sametet
        referencestring = 'sametet';
    elseif statictet
        referencestring = sprintf('staticref @ %s',statictet_region);       % static reference
    end

end
if runscript
    animals_torun =  {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
end
if postprocessing
    animals_toprocess =  {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Module 0: Plot a single selected units' theta phase histograms

if plot_results
    
    % unit specification
    animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
    %manual_adtc = {sortrows(classtable6_adtc{2}{1},[1 2 3 4]),[1 .5 .5],'CA2 P'};    % P
    %manual_adtc = {sortrows(classtable6_adtc{1}{1},[1 2 3 4]),[.2 .2 .2],'CA1'};    % CA1
    %manual_adtc = {sortrows(classtable6_adtc{3}{1},[1 2 3 4]),[1 .3 .3],'CA3'};    % CA3
    load('/opt/data13/kkay/Superlin_data/Classtable_30-Dec-2014.mat')
    load('/opt/data13/kkay/Superwell_data/wellunitlist_01-Jan-2015.mat')
    
    %manual_adtc = {sortrows([classtable6_adtc{2}{2} ; classtable6_adtc{2}{3} ],[1 2 3 4]),[.6 .7 1],'CA2 N'};    % CA2 N
    manual_adtc = {sortrows([classtable6_adtc{2}{1} ; classtable6_adtc{2}{4} ],[1 2 3 4]),[1 .5 .5],'CA2 P'};    % CA2 P
    %manual_adtc = {sortrows([classtable6_adtc{3}{1} ; classtable6_adtc{3}{2}],[1 2 3 4]),[1 .3 .3],'CA3'};    % CA3
    %manual_adtc = {sortrows([classtable6_adtc{1}{1} ; classtable6_adtc{1}{2}],[1 2 3 4]),[.2 .2 .2],'CA1'};    % CA1
    
    %manual_adtc = {sortrows([ca1_well],[1 2 3 4]),[.2 .2 .2],'CA1'};    % CA1 well
    %animals_toplot = {'Government','Egypt','Chapati','Higgs','Frank'} ;  % filter manual_adtc for these animals only
    animals_toplot = {'Government','Egypt','Chapati','Dave','Higgs','Frank'}; 
    %animals_toplot = {'Dave'} ;  % filter manual_adtc for these animals only
        
    if plot_singles
        min_spikes = 20;
        normhist = 1;
        maxunitsinfig = 10;
        manual_clr = manual_adtc{2};
    end
     
    if plot_population
        min_numspikes = 50;
        % plot row 4 : moddepth histc edges
        moddepth_edges = 0:.01:1;
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if runscript
    
    clear f
 
    iterator = 'singlecelleeganal';
    
    for animal = animals_torun
        
        % reference tetrode selection
        if sametet && ~statictet
            eegfilter = {'geteegtet', 'theta', 'sametet',1};
        elseif statictet && ~sametet 
            statictet_tetnum = nan;
            if strcmp(statictet_region,'CA1')
                switch animal{1}
                    case 'Conley'
                        statictet_tetnum = 30;
                    case 'Higgs'
                        statictet_tetnum = 21;  % 21 is a CA1b tetrode that gives CA3 trough and CA1 (few units) trough
                    case 'Government'
                        statictet_tetnum = 17;  % callosum tetrode:  matches Csicsvari for CA1 firing phase
                    case 'Chapati'
                        statictet_tetnum = 1;   % callosum tetrode:  matches Csicsvari for CA1 firing phase
                    case 'Egypt'
                        statictet_tetnum = 11;  % callosum tetrode: matches Csicsvari for CA1 firing phase
                    case 'Frank'
                        statictet_tetnum = 30;  % "Reference" (callosum) tetrode: matches Csicsvari for CA1 firing OK -- could use a pyramidal that has better agreement
                    case 'Bond'
                        statictet_tetnum = 30;
                    case 'Corriander'
                        statictet_tetnum = 24;
                    case 'Apollo'
                        statictet_tetnum = 8;      % corpus callosum
                    case 'Dave'
                        statictet_tetnum = 1;      % corpus calllosum
                        
                        %                 case 'Barack'
                        %                     statictet_tetnum = 4; % Barack: reference is 4 for days 1-17
                        %                 case 'Calvin'
                        %                     statictet_tetnum = 16; % Calvin: reference is 16
                        %                 case 'Dwight'
                        %                     statictet_tetnum = 7; % Dwight; reference is 7
                end
            elseif strcmp(statictet_region,'CA3b')
                switch animal{1}
                    %case 'Higgs'
                    %    statictet_tetnum = 12;   % CA3c <<< with lots of units
                    case 'Government'
                        statictet_tetnum = 13;   % consistently the most units
                    case 'Chapati'
                        statictet_tetnum = 4;    % best guess for now
                    case 'Egypt'
                        statictet_tetnum = 14;   % consistently the most units throughout recording
                    case 'Frank'
                        statictet_tetnum = 29;   % 7, 8, 29 look depth-constant + have of units
                end
            end 
            eegfilter = {'geteegtet', 'theta', 'statictet',statictet_tetnum};
        else
            error('must choose either sametet or statictet')
        end        
        
        
        disp(['*************** ANIMAL ' animal{1} ' *************************'])
        % grab cells inclusively, then sort their identities later in postprocessing
        % note that for animals in which .hemisphere field exists,
        % get cells just from the Right hippocampus
        clear cellfilter
        if strcmp(animal{1},'Frank') || strcmp(animal{1},'Bond')
            cellfilter = '(isequal($type, ''principal'')) && (~isequal($hemisphere, ''left''))';
        else
            cellfilter = '(isequal($type, ''principal''))';
        end
        
        for state = 1:length(statespec)
          
            timefilter = statespec{state}{1};
            epochfilter = epochmaker(statespec{state}{2});
            
            f{1}{state} = createfilter('animal',animal{1},'epochs',epochfilter, 'excludetimefilter', timefilter, 'cells',cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
            
            if statictet
                f{1}{state} = setfilterfunction(f{1}{state}, 'dfakk_thetaphase', {'spikes', 'theta'},'nbins',24);
            elseif sametet
                f{1}{state} = setfilterfunction(f{1}{state}, 'dfakk_thetaphase', {'spikes', 'theta'},'nbins',24);
                %f{1}{state} = setfilterfunction(f{1}{state}, 'dfakk_thetagndphase', {'spikes', 'thetagnd'},'nbins',24);
            end
            
            f{1}{state} = runfilter(f{1}{state});
            
            % append runscript parameters to the raw output structs
            datestring = date;
            runscript_params = paramsstruct(timefilter,epochfilter,datestring);
            f{1}{state}.runscript_params = runscript_params;            
        
        end
        disp(['done w/ ' animal{1}] )
        save(sprintf('%sSupertheta_raw_%s',datadir,animal{1}(1:3)),'f','-v7.3');
        clear f;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Consolidate single cells' across epochs in a day (.celloutput field)

if postprocessing
    
    clear supertheta
    
    for animal = animals_toprocess
        
        supertheta.animalname = animal{1};
        supertheta.bins = bins;
        supertheta.dtc = [];
        supertheta.data = [];       
        
        % load data
        filename = dir(sprintf('%sSupertheta_raw_%s*',datadir,animal{1}(1:3)));
        load([datadir filename.name])        
        
        % collect detc indices from state 1 
        detc = [];
        for kk=1:length(f{1}{1}.output{1})
            detc = [detc; f{1}{1}.output{1}(kk).index];
        end
        dtc = unique(detc(:,[1 3 4]),'rows');
        supertheta.dtc = dtc;
        
        % iterate through epoch data blocks
        clear C;
        
        for kk=1:size(dtc,1)
            
            % initialize struct output
            C.dtc = dtc(kk,:);
            C.epochs = [];
            C.spikephases = {};
            C.numspikes = zeros(1,length(statespec));
            C.included_duration = zeros(1,length(statespec));
            C.meanphase = nan(1,length(statespec));
            C.moddepth = nan(1,length(statespec));
            C.phasehist = cell(1,length(statespec));
            C.rayleighp = nan(1,length(statespec));
            
            % first find epochs that match cell
            ind = [];
            epochs = [];
            while rowfind(dtc(kk,:),detc(:,[1 3 4])) ~= 0
                row = rowfind(dtc(kk,:),detc(:,[1 3 4]));
                ind = [ind row];
                epochs = [ epochs  detc(row,2)  ];
                detc(row,:)=[nan nan nan nan];
            end
            C.epochs = epochs;
            
            for state = 1:length(statespec)
                % retrieve phase data from each epoch
                C.spikephases{state} = [];
                for rr=ind
                    dataentry = f{1}{state}.output{1}(rr);
                    C.spikephases{state} = [C.spikephases{state} ; dataentry.sph];
                    C.numspikes(state) = C.numspikes(state) + length(dataentry.sph);
                    C.included_duration(state) = C.included_duration(state) + dataentry.included_duration;
                end
                if ~isempty(C.spikephases{state})
                    % compute circular stats
                    clear i;
                    A=mean(exp(i*C.spikephases{state}));
                        meancos=real(A);
                        meansin=imag(A);
                    C.meanphase(state)=atan2(meansin,meancos);
                    C.moddepth(state)=abs(A);
                    C.phasehist{state}=histc(C.spikephases{state},bins);
                    % circular significance
                    out=rayleigh_test(C.spikephases{state});
                    C.rayleighp(state)=out.p;
                else
                    C.meanphase(state)=nan;
                    C.moddepth(state)=nan;
                    C.phasehist{state}=nan;
                    C.rayleighp(state)=nan;
                end
            end
            
            supertheta.data = [supertheta.data   C]  ;
        end
        
        save(sprintf('%sSupertheta_%s',datadir,animal{1}(1:3)),'supertheta','-v7.3');
        
        clear supertheta; clear f;
    end
end

%% Module 0: Plot histogram of manually specified units

if plot_results && plot_singles
    
    % first filter manual_adtc for animals to plot %%%%%%%%%%%%%%%%%
    anim = []; animstring = '';
    for aa = 1:length(animals_toplot)
    	anim = [anim find(strcmp(animals_toplot{aa},animals_order))];
        animstring = [animstring  animals_toplot{aa}(1:3)];
    end
    manual_adtc{1} = manual_adtc{1}(ismember(manual_adtc{1}(:,1),anim),:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    
    % first sort manual_adtc
    adtc_toplot = sortrows(manual_adtc{1},[1 2 3 4])
    % next determine which animals to load
    animals_toload = animals_order(unique(adtc_toplot(:,1))')

    unitcount = 0;  % # of cells in figure window
    K = figure('units','normalized','outerposition',[0 0 .5 1]);
    
    for animal = animals_toload
        
        animalnum = find(strcmp(animal,animals_order));
        dtc_animal = adtc_toplot(adtc_toplot(:,1)==animalnum,2:4);
        % load data
        filename = dir(sprintf('%sSupertheta_%s*',datadir,animal{1}(1:3)));
        load([datadir filename.name])
        for cc = 1:size(dtc_animal,1)
            dtc = dtc_animal(cc,:);
            % get unit's data
            ind = rowfind(dtc,supertheta.dtc);
            if ~ind
                disp(sprintf('%s %s not found in supertheta',animal,num2str(dtc)));
            else
                dataentry = supertheta.data(ind);
                numstates = length(dataentry.numspikes);
                unitcount = unitcount + 1;
                if unitcount > maxunitsinfig
                    K = figure('units','normalized','outerposition',[0 0 .5 1]);
                    unitcount = 1;
                end  
            end
            for state = 1:numstates
                numspikes = dataentry.numspikes(state);
                subplot(maxunitsinfig,numstates,(unitcount-1)*numstates + state)
                % plot
                N = histc(dataentry.spikephases{state},supertheta.bins);
                    N=N(1:(end-1));
                    N=N(:);
                if normhist
                    N = N/sum(N);                    
                end
                bins_plot = bins(1:(end-1));
                bins_plot = bins_plot + (bins(2)-bins(1))/2;
                h = bar([bins_plot bins_plot+2*pi],[N ; N],'histc');
                set(h(1),'facecolor',manual_clr); set(h(1),'edgecolor','none');
                if state == 1
                    title(sprintf('%s %s (%d)',animal{1}(1:3),num2str(dtc),numspikes),'FontSize',14,'FontWeight','bold')
                else
                    title(sprintf('(%d)',numspikes),'FontSize',14,'FontWeight','bold')                    
                end
                axis tight
                hold on
                % guide lines
                plot([pi,pi],[0 2*max(N)],'k--','LineWidth',1.5)
                plot([-pi,-pi],[0 2*max(N)],'k--','LineWidth',1.5)
                plot([3*pi,3*pi],[0 2*max(N)],'k--','LineWidth',1.5)
                plot([0,0],[0 2*max(N)],'k:','LineWidth',1.5)
                plot([2*pi,2*pi],[0 2*max(N)],'k:','LineWidth',1.5)
                clear ylim;
                if ~normhist
                    if max(N) > 0
                        ylim([0 1.5*max(N)])
                    end
                else
                    ylim([0 .1])
                end
                statecolor
                set(gca,'color',color)
            end
        end
    end
end

%% Module 0: Plot phase histogram of manually specified units

if plot_results && plot_population
    
    groupname = manual_adtc{3};
    groupcolor = manual_adtc{2};
    
    % first filter manual_adtc for animals to plot %%%%%%%%%%%%%%%%%
    anim = []; animstring = '';
    for aa = 1:length(animals_toplot)
    	anim = [anim find(strcmp(animals_toplot{aa},animals_order))];
        animstring = [animstring  animals_toplot{aa}(1:3)];
    end
    manual_adtc{1} = manual_adtc{1}(ismember(manual_adtc{1}(:,1),anim),:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % sort
    adtc_toplot = sortrows(manual_adtc{1},[1 2 3 4]);
    % determine which animals to load
    animals_toload = animals_order(unique(adtc_toplot(:,1))');

    % data to collect
    unitcount = zeros(1,numstates);     % # of cells with minimum spikes
    sigcount = zeros(1,numstates);
    adtc_minspike = cell(1,numstates);  % indices of units with minimum spikes
    adtc_sigmod = cell(1,numstates);    % indices of units that have rayleigh < 0.05
    histdata_raw = cell(1,numstates);       % storage of raw psths
    histdata_norm = cell(1,numstates);
    histdata_meanphases = cell(1,numstates);
    MODDEPTHS = cell(1,numstates);      % storage of moddepths
    
    % first collect data
    for animal = animals_toload
        
        animalnum = find(strcmp(animal,animals_order));
        dtc_animal = adtc_toplot(adtc_toplot(:,1)==animalnum,2:4);
        % load data
        filename = dir(sprintf('%sSupertheta_%s*',datadir,animal{1}(1:3)));
        load([datadir filename.name])
        for cc = 1:size(dtc_animal,1)
            dtc = dtc_animal(cc,:);
            % get unit's data
            ind = rowfind(dtc,supertheta.dtc);
            if ~ind
                disp(sprintf('%s %s not found in supertheta',animal,num2str(dtc)));
            else
                dataentry = supertheta.data(ind);
                numstates = length(dataentry.numspikes);
            end
            for state = 1:numstates
                
                if dataentry.numspikes(state) < min_numspikes
                    continue
                else
                    adtc_minspike{state} = [adtc_minspike{state} ; animalnum dtc];
                    unitcount(state) = unitcount(state) + 1;
                    MODDEPTHS{state} = [MODDEPTHS{state} ; dataentry.moddepth(state)];   % taken for all units, regardless of significance
                    % calculate histogram  (for all units)
                    N = histc(dataentry.spikephases{state},supertheta.bins);
                        N=N(1:(end-1));
                        N=N(:)';
                    histdata_raw{state} = [histdata_raw{state} ; N];
                    histdata_norm{state} = [histdata_norm{state} ; N/sum(N)];
                    % if significant, get mean phase
                    if dataentry.rayleighp(state) < 0.05
                        adtc_sigmod{state} = [adtc_sigmod{state} ; animalnum dtc];
                        sigcount(state) = sigcount(state) + 1;
                        histdata_meanphases{state} = [histdata_meanphases{state} ; dataentry.meanphase(state)];
                    end                    
                end
            end
        end
    end
    
    
    figure('units','normalized','outerposition',[0 0 .5 1]);

    % plot #1: normhist
    for state = 1:numstates

        
        subplot(4,numstates,state)
        
        bins_plot = bins(1:(end-1));
        bins_plot = bins_plot + (bins(2)-bins(1))/2;
        % if normhist, calculate mean and SEM hist
            meanhist = mean(histdata_norm{state},1)';
            semhist = std(histdata_norm{state},[],1)/sqrt(size(histdata_norm{state},1));
            numunits = size(histdata_norm{state},1);
            if isempty(meanhist)
                continue
            end
            % plot hists
            h = bar([bins_plot bins_plot+2*pi],[meanhist ; meanhist],'histc');
            set(h(1),'facecolor',groupcolor); set(h(1),'edgecolor','none');
            hold on
            % plot sem 
           for jj=1:length(bins_plot)
                plot([bins_plot(jj),bins_plot(jj)],[meanhist(jj)-semhist(jj) meanhist(jj)+semhist(jj)],'k','LineWidth',2)
                plot([bins_plot(jj)+2*pi,bins_plot(jj)+2*pi],[meanhist(jj)-semhist(jj) meanhist(jj)+semhist(jj)],'k','LineWidth',2)
            end            
            clear ylim; axis tight;
            maxval = 0.06;
            ylim([0 maxval])
           
        if state == 1
            title(sprintf('n = %d',numunits),'FontSize',14,'FontWeight','bold')
        else
            title(sprintf('n = %d',numunits),'FontSize',14,'FontWeight','bold')
        end
        hold on
        % guide lines
        plot([pi,pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([-pi,-pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([3*pi,3*pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([0,0],[0 2*maxval],'k:','LineWidth',1.5)
        plot([2*pi,2*pi],[0 2*maxval],'k:','LineWidth',1.5)
        
        statecolor
        set(gca,'color',color)

    end
    
    
    % plot #2 aggregate spike hist
    for state = 1:numstates

        subplot(4,numstates,numstates + state)
        
        statecolor
        
        bins_plot = bins(1:(end-1));
        bins_plot = bins_plot + (bins(2)-bins(1))/2;
        % if normhist, calculate mean and SEM hist
            numunits = size(histdata_raw{state},1);
            sumhist = sum(histdata_raw{state},1)';
            h = bar([bins_plot bins_plot+2*pi],[sumhist ; sumhist],'histc');
            set(h(1),'facecolor',groupcolor); set(h(1),'edgecolor','none');
            clear ylim; axis tight;
            ylim([0 1.5*max(sumhist)])
            maxval = max(sumhist);
    
            title(sprintf('n = %d',numunits),'FontSize',14,'FontWeight','bold')

        hold on
        % guide lines
        plot([pi,pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([-pi,-pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([3*pi,3*pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([0,0],[0 2*maxval],'k:','LineWidth',1.5)
        plot([2*pi,2*pi],[0 2*maxval],'k:','LineWidth',1.5)
                statecolor
        set(gca,'color',color)
    end    

    % plot #3 mean phase of units
    for state = 1:numstates

        subplot(4,numstates,numstates*2 + state)
        
        statecolor
        
        bins_plot = bins(1:(end-1));
        bins_plot = bins_plot + (bins(2)-bins(1))/2;
        % calculate hist over all units & plot
        N = histc(histdata_meanphases{state},supertheta.bins);
            N=N(1:(end-1));
            N=N(:)';       
            if isempty(N)
                continue
            end
            h = bar([bins_plot bins_plot+2*pi],[N  N],'histc');
            set(h(1),'facecolor',groupcolor); set(h(1),'edgecolor','none');
            clear ylim; axis tight;
            ylim([0 1.5*max(N)])
            maxval = max(N);
    
        % title
        title(sprintf('n = %d',sigcount(state)),'FontSize',14,'FontWeight','bold')

        hold on
        % guide lines
        plot([pi,pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([-pi,-pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([3*pi,3*pi],[0 2*maxval],'k--','LineWidth',1.5)
        plot([0,0],[0 2*maxval],'k:','LineWidth',1.5)
        plot([2*pi,2*pi],[0 2*maxval],'k:','LineWidth',1.5)
                statecolor
        set(gca,'color',color)
    end        
    
    
    
    % plot #4 moddepths
    for state = 1:numstates
        
        subplot(4,numstates, numstates*3 + state)
        
        statecolor
        
        N = histc(MODDEPTHS{state},moddepth_edges);
        if isempty(N)
            continue
        end
        h = bar(moddepth_edges,N,'histc');
        xlim([0 .7]);
        set(h(1),'facecolor',groupcolor); set(h(1),'edgecolor','none');   
        
        title(sprintf('n = %d',length(MODDEPTHS{state})),'FontSize',14,'FontWeight','bold')
                statecolor
        set(gca,'color',color)
    end
    
    % supertitle
    titlestring = sprintf('%s theta (%s, ref: %s)',groupname,animstring,referencestring);
        [~,title_handle] = suplabel(titlestring,'t');
        set(title_handle,'Fontsize',18,'FontWeight','bold','Color','k')
        %set(gcf, 'renderer', 'zbuffer')     
    
    
end



