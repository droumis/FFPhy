datadir = '/opt/data13/kkay/Superclustdecode_data/100ms_decodes/';



setbands = 1;
calculate = 0;   % loads P file and identify candidate non-SWR bins, updates P file
powercalc = 1;  % loads P files and for each bin retrieves
                 % difference in distance between decoded and actual +
                 % evaluates gamma power + ripple power
basicstudy = 1;

    if setbands
        clear tracespec
        tracespec{1} = {'rippletrace',1,[.5 .5 .5]};   % tracefile name + TFgroup
        tracespec{2} = {'lowgammatrace',1,[1 0 0]};
        tracespec{3} = {'fastgammatrace',2,[.4 .6 1]};
    end

    if calculate 
       animals_tocalc_1 = {'Bond','Frank','Government','Dave','Corriander'};  %'Egypt','Government','Chapati','Higgs',
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       manual_epochs = [];    
       epochfilter = epochmaker('runW');
       % target times
       timefilterscript
       timefilter_target = { immobile2 };
       % remove bins that overlap with 2 SD ripples
       remove_rip = 1;
       % maximum ripple power
       maxrippower = 0.25;  % remove bins that have +0.5 SD or more ripple power
       save_update_P = 1;  % updates P file with the nonlocal bins in .nonlocalbins
    end
      
    
    if powercalc 
       animals_tocalc_2 =  {'Dave'};  %{'Frank'}; %{'Government','Dave','Bond','Corriander'};  % ,'Frank'  ,'Frank'
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       manual_epochs = [];    
       ignore_days = [];
       epochfilter = epochmaker('runW');     
       kern_ms = [];
    end
               
                 
  if basicstudy
        approach1 = 1;      % local vs. nonlocal, use a location threshold
            locthresh = 2;
            boxstyle = 0;
        approach2 = 0;      % WG vs. LG
        approach3 = 0;      % correlate diffdist w/ LG-WG, Spearman's corr
  end
  

        
    
    
    
    if calculate
        
        for aa = 1:length(animals_tocalc_1)
            animalname = animals_tocalc_1{aa};
            animalinfo = animaldef(animalname);
            
            animpref = animalname(1:3);
            animalinfo = animaldef(animalname);
            daydir = getdaydir(animalname);
            
            task = loaddatastruct(animalinfo{2},animalinfo{3},'task');
            
            if ~isempty(manual_epochs)
                epochs = manual_epochs;
            else
                epochs = evaluatefilter(task,epochfilter);
            end
            
            for ee = 1:size(epochs,1)
                
                d = epochs(ee,1);
                ep = epochs(ee,2);
                disp([ animpref ' : ' num2str([d ep])])
                
                pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',d);
                postimevec = pos{d}{ep}.data(:,1);
                
                % Identify target periods for detecting candidate bins
                [targout,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},timefilter_target,[d ep]);
                targetperiods = targout{d}{ep};
                
                % Load decode data
                if ~exist('P','var') || ~strcmp(P.animalname,animalname) || ...
                        ~all(P.dayep == [d ep])
                    cd(datadir)
                    load(sprintf('%s_%d_%d.mat',animpref,d,ep),'P');
                end
                
                % Basic values
                numbins = size(P.posteriors{1},2);
                activebins = P.activebins{1};
                if isfield(P,'centerarmmax')
                    centerarmmax = P.centerarmmax(1);
                    rightarmmax = P.rightarmmax(1);
                    binvec_c = P.binvec_c{1};
                    times_activebins = P.binvec_c{1}(activebins);
                else
                    centerarmmax = P.lastcenterbin(1);
                    rightarmmax = P.lastrightbin(1);     
                    binvec_c = P.binvec_c;
                    times_activebins = P.binvec_c(activebins);
                end

                binsize = binvec_c(2) - binvec_c(1);
                
                % Initialize outputs
                armdists = cell(1,2);  % 1: actual animal posn
                                       % 2: decoded posn: peak of posterior
                                       % 3: decoded posn: mean of posterior

                % Obtain position of animal in each bin (from bincenter time)
                armdists{1} = P.armdists(lookup(binvec_c,postimevec));
                % peak of posteriors
                [~,maxind] = max(P.posteriors{1},[],1);
                armdists{2} = P.linposbins{1}(maxind);
%                 % mean of posterior  *** needs circular treatment at center junction, not completed here..
%                 postr = P.posteriors{1};
%                 postsums = sum(postr,1);
%                 decpos = nan(1,length(postsums));  % 1st moment decoded position
%                 for ll = 1:length(postsums)
%                     postr(:,ll) = postr(:,ll)/postsums(ll) ;
%                     decpos(ll) = sum(postr(:,ll) .* P.linposbins{1});
%                 end
%                 armdists{3} = decpos;

                
                % For each position condition, process arm position
                armnum = cell(1,2);
                armpos = cell(1,2);
                for xx = 1:2
                    armnum{xx} = zeros(1,numbins);
                    armpos{xx} = zeros(1,numbins);
                    % identify which arm the location is on
                    centerind = armdists{xx} <= centerarmmax;
                    rightind = (armdists{xx} > centerarmmax) & (armdists{xx} <= (centerarmmax + rightarmmax));
                    leftind = armdists{xx} > (centerarmmax + rightarmmax);
                    armnum{xx}(centerind) = 1;  % center
                    armnum{xx}(rightind) = 2;  % right
                    armnum{xx}(leftind) = 3;  % left
                    % obtain distance along each respective arm ( *** with the center junction as the 0, "origin")
                    armpos{xx}(centerind) = centerarmmax - armdists{xx}(centerind);             % center
                    armpos{xx}(rightind) = armdists{xx}(rightind) - centerarmmax;               % right
                    armpos{xx}(leftind) = armdists{xx}(leftind) - centerarmmax - rightarmmax;   % left
                    % sanity check
                    if 0
                       indsums = sum(centerind) + sum(rightind) + sum(leftind);
                       disp(sprintf('%d == %d inds',length(armdists{xx}),indsums)) 
                       if any(armpos{xx} < 0)
                           keyboard
                       end
                    end
                end
                
                % Calculate distance difference
                diffdist = nan(1,numbins);
                sepinds = armnum{1} ~= armnum{2};   % separate arms
                sameinds = armnum{1} == armnum{2};  % same arms
                diffdist(sepinds) = armpos{1}(sepinds) + armpos{2}(sepinds);  % if separate, then add respective arm positions
                diffdist(sameinds) = abs( armpos{1}(sameinds) - armpos{2}(sameinds) );
                
                if 0
                    targtime = 784.5;
                    targbin = lookup(targtime,binvec_c);
                    diffdist(targbin)
                    keyboard
                end
                
                candinds = [];
                numkilled = 0;
                
                % First filter for bins that fall within target periods
                
                    targetinds = logical(isExcluded(times_activebins,targetperiods))';
                        candinds = activebins(targetinds);
                    
                % Second, remove bins that contain or overlap with 2 SD SWRs
                if remove_rip
                    % Load SWR data
                    timefilterscript
                    outperiods = evaluatetimefilter(animalinfo{2},animalinfo{3},{rip2},[d ep]);
                    ripperiods = outperiods{d}{ep};
                    for nn = length(candinds):-1:1
                        j = candinds(nn);
                        bin_start = binvec_c(j) - binsize/2;
                        bin_end = binvec_c(j) + binsize/2;
                        ripoverlap =  logical(  isExcluded(bin_start,ripperiods) || isExcluded(bin_end,ripperiods) || ...
                                                any(isExcluded(ripperiods(:),[bin_start bin_end]))   ) ;
                        if ripoverlap
                            candinds(nn) = [];
                            numkilled = numkilled + 1;
                            %disp('rip removal')
                        end
                    end
                end                        
                        
                % Third, if specified, filter for bins that don't exceed given max ripple power        
                if ~isempty(maxrippower)
                    % load swr data first
                    riptrace = loadtracestruct(animalinfo{2},animalinfo{3},'rippletrace',d,ep);
                    tvec = riptrace{d}{ep}{1}.eegtimesvec_ref;
                    mean_runimmo = riptrace{d}{ep}{1}.mean_runimmo;
                    stdev_runimmo = riptrace{d}{ep}{1}.stdev_runimmo;
                    ztrace = (riptrace{d}{ep}{1}.powertrace - mean_runimmo) / stdev_runimmo ;
                    meanpow = nan(1,length(candinds));
                    for bb = 1:length(candinds)
                         k = candinds(bb);
                         if k-1 <= 0 || k + 1 >= length(binvec_c)
                             meanpow(bb) = 9999;
                             break
                         end
                         start_bin = mean([binvec_c(k-1) binvec_c(k)]);
                         end_bin = mean([binvec_c(k) binvec_c(k+1)]);
                         f = lookup(start_bin,tvec);
                         g = lookup(end_bin,tvec);
                         meanpow(bb) = mean(ztrace(f:g));
                    end
                    killbins = meanpow > maxrippower;
                    % remove any bins that are contained
                    numkilled = numkilled + sum(killbins);
                    candinds(killbins) = [];
                end
                
                
                disp(sprintf('%d bins detected (%d above rip power)',length(candinds),numkilled))
                
                if save_update_P
                    P.candbins = candinds;
                    P.diffdist = diffdist;   % difference in distance between actual and decoded
                    % save output
                    cd(datadir)
                    disp('update-saving non local bins in P file')
                    savefilename = sprintf('%s_%d_%d',animalname(1:3),d,ep);
                    save(savefilename,'P')
                    clear nonlocalbins
                end
            end
        end
    end
    
    
    
    
    
    if powercalc
        
        bindata = {};  % {animalnum}[ <diff dist>  <ripplez>  <lowgammaz>  <fastgammaz>   ];
        numbands = length(tracespec);
        
        for aa = 1:length(animals_tocalc_2)
            
            bindata{aa} = [];
            
            animalname = animals_tocalc_2{aa};
            animalinfo = animaldef(animalname);
            
            animpref = animalname(1:3);
            animalinfo = animaldef(animalname);
            daydir = getdaydir(animalname);
            
            task = loaddatastruct(animalinfo{2},animalinfo{3},'task');
            
            if ~isempty(manual_epochs)
                epochs = manual_epochs;
            else
                epochs = evaluatefilter(task,epochfilter);
            end
            
            for ee = 1:size(epochs,1)
                
                d = epochs(ee,1);
                ep = epochs(ee,2);
                
                if ~isempty(ignore_days)
                    if ismember(d,ignore_days)
                        disp('day ignored')
                        continue
                    end
                end
               
                % Load decode data
                cd(datadir)
                load(sprintf('%s_%d_%d.mat',animpref,d,ep),'P');
                
                numcandbins = length(P.candbins);
                
                if numcandbins > 0
                    
                    % Get distance difference (decoded from actual position)
                    diffdist = P.diffdist(P.candbins)';
                    
                    % Now obtain gamma and ripple powers (from 4 ms smoothed trace data)
                    binpowers = nan(numcandbins,numbands);
                    
                    for g = 1:numbands
                        
                        % Load power trace data
                        out = loadtracestruct(animalinfo{2}, animalinfo{3},tracespec{g}{1}, d, ep);
                        TF = tracespec{g}{2};
                        if 1
                            zpowertrace = zscoretrace(out{d}{ep}{TF},kern_ms);
                        elseif 1
                            kernsamps = 1500 * (kern_ms)/1000;
                            gauskern = gaussian(kernsamps,8*kernsamps);
                            zpowertrace = zscorer(smoothvect(out{d}{ep}{TF}.powertrace,gauskern));
                        end
                        timevec = out{d}{ep}{TF}.eegtimesvec_ref;
                        
                        for dd = 1:numcandbins
                            k = P.candbins(dd);
                            if isfield(P,'centerarmmax')
                                start_bin = mean([P.binvec_c{1}(k-1) P.binvec_c{1}(k)]);
                                end_bin = mean([P.binvec_c{1}(k) P.binvec_c{1}(k+1)]);
                            else
                                start_bin = mean([P.binvec_c(k-1) P.binvec_c(k)]);
                                end_bin = mean([P.binvec_c(k) P.binvec_c(k+1)]);                                
                            end
                            a = lookup(start_bin,timevec);
                            b = lookup(end_bin,timevec);
                            binpowers(dd,g) = mean(zpowertrace(a:b));
                        end
                        
                    end
                    
                    bindata{aa} = [bindata{aa} ; diffdist binpowers];
                    
                end
                
            end
            
            totalcandbins = size(bindata{aa},1);
            
            disp(sprintf('%s : %d total cand bins',animalinfo{3},totalcandbins))
            
        end
        
        
        
        
    end
    
    

if basicstudy
    
    alldata = [];
    for an = 1:length(bindata)
        alldata = [alldata ; bindata{an}];
    end
    
    % Histogram all decoded bins to check for nonlocality
    if 1
        figure;
        hist(alldata(:,1),100);
        titlestring = sprintf('Decode diff, all candidate bins (n = %d, %d animals)',...
            size(alldata,1),length(bindata));
        title(titlestring,'fontsize',14,'fontweight','bold');
    end
    
    % Approach 1: classify bins into local vs. nonlocal
    if approach1
       
        H = figure('units','normalized','outerposition',[0 0 .8 .3]);
        
        clear inds
        inds{1} = alldata(:,1) < locthresh;
        inds{2} = alldata(:,1) > locthresh;
            numlocbins = sum(inds{1});
            numnonlocbins = sum(inds{2});
            
        % direct band power comparison
        for b = 1:3
            subplot(1,4,b)
            clr = tracespec{b}{3};
            
            if ~boxstyle
                for l = 1:2
                    meanpow = mean(alldata(inds{l},b+1));
                    sempow = std(alldata(inds{l},b+1))/sqrt(sum(inds{l}));
                    h = bar(l,meanpow); hold on ; set(h,'BarWidth',2,'facecolor',clr);
                    errorbar(l, meanpow, sempow, 'k', 'linestyle', 'none','linewidth',2);
                end
            else
                vals1 = alldata(inds{1},b+1) ;
                vals2 = alldata(inds{2},b+1) ;
                boxplot([vals1; vals2],[ones(length(vals1),1) ; 2 * ones(length(vals2),1)],'notch','on');
            end
            
            set(gca,'fontsize',12);
            if b == 1
            title(sprintf(' %s pow (n = %d, %d)',...
                          tracespec{b}{1}(1:end-5),numlocbins,numnonlocbins),'fontweight','bold','fontsize',14)
            else
            title(sprintf(' %s pow',...
                          tracespec{b}{1}(1:end-5)),'fontweight','bold','fontsize',14)                
            end
            set(gca,'xtick',1:2,'XTickLabel',{'local','nonlocal'})          
            ylim([-.2 .1])
        end
        
        % difference in lowgamma vs. fastgamma
        
         subplot(1,4,4)
        if boxstyle
            vals1 = alldata(inds{1},3) - alldata(inds{1},4);
            vals2 = alldata(inds{2},3) - alldata(inds{2},4);
            boxplot([vals1; vals2],[ones(length(vals1),1) ; 2 * ones(length(vals2),1)],'notch','on');
        else
            for l = 1:2
                vals = alldata(inds{l},3) - alldata(inds{l},4);
                meandiff = nanmean(vals);
                semdiff = nanstd(vals)/sqrt(length(vals));
                h = bar(l,meandiff); hold on ; set(h,'BarWidth',2,'facecolor',[.5 1 .4]);
                errorbar(l, meandiff, semdiff, 'k', 'linestyle', 'none','linewidth',2);
            end
        end
        set(gca,'XTickLabel',{'local','nonlocal'})
        set(gca,'fontsize',12);
        title(sprintf('LG - FG diff (n = %d, %d)',numlocbins,numnonlocbins),'fontweight','bold','fontsize',14)
        set(gca,'xtick',1:2,'XTickLabel',{'local','nonlocal'})   
        ylim([-.2 .1])
        
    end
    
    % Approach 2: classify bins as either WG or LG dominated
    if approach2
        
        clear inds
        gammadiff = alldata(:,3) - alldata(:,4) ;
        inds{1} =  gammadiff > 0;  % LG bins
        inds{2} = gammadiff < 0;  % WG bins
        
        figure
        
        for ff = 1:2
            if ff == 1
                clr = [1 0 0];
            elseif ff == 2
                clr = [.4 .6 1];
            end
            meandist = mean(alldata(inds{ff},1));
            semdist = std(alldata(inds{ff},1))/sqrt(sum(inds{ff}));
            h = bar(ff,meandist); hold on ;
            set(h,'BarWidth',2,'facecolor',clr);
            errorbar(ff, meandist, semdist, 'k', 'linestyle', 'none','linewidth',2);
            set(gca,'XTickLabel',{'LG','FG'})
            set(gca,'fontsize',12);
            title(sprintf('Decoded Excursion: LG vs. FG'),'fontweight','bold','fontsize',16)
        end
        
       
    end    
    
    % Approach 3: Correlations
    if approach3
        
        gammadiff = alldata(:,3) - alldata(:,4) ;
        
        figure
        %scatter(alldata(:,1),gammadiff,20,'k','o')
        scatter(gammadiff,alldata(:,1),20,'k','o')
        [RHO,PVAL] = corr(alldata(:,1),alldata(:,4),'type','Spearman')
    end      
    
    
end