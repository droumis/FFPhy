% studies the timescale at which two different gammas correlate

calculate = 0;
    if calculate
        plot_test = 0;
    end
plot_animal = 1;
    if plot_animal
        if 1
            animals_toplot = {'Government'}; %s
        else
            animals_toplot = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
        end
        datadir = '/opt/data13/kkay/__WG/Super_running_corr/Runepochs';
        plot_xlimit = [0 5000];
        plot_allepochs = 1;
    end
    
if calculate
    animals_tocalc =  {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
    %animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
    corrwin = 0.5;  % in sec
    kernsigs = [10:10:300 350:50:1000 1250:250:5000];  % in ms
        numkernsigs = length(kernsigs);
    clear tracespec;
    tracespec{1} = {'fastgammatrace',2,[.4 .6 1]};   % tracefile name + TFgroup
    tracespec{2} = {'lowgammatrace',1,[1 0 0 ]};
    timefilterscript;   clear statespec;
    if 0
        % RUN
        epochstring = 'runW';
        statespec{1} = { fullepoch } ;
        statespec{2} = { vel10 } ;
        statespec{3} = { immobile2  } ;
        statespec{4} = { vel4  } ;
        datadir = '/opt/data13/kkay/__WG/Super_running_corr/Runepochs';
    elseif 1
        % SLEEP
        epochstring = 'sleep';
        statespec{1} = { waking } ;
        statespec{2} = { sleepc, noREM } ;
        datadir = '/opt/data13/kkay/__WG/Super_running_corr/Sleep_wake_sleepc';
    end
    numstates = length(statespec);
end


if calculate
    
    A = [];
    
    for aa = 1:length(animals_tocalc)
        
        animalname = animals_tocalc{aa};
        animalinfo = animaldef(animalname);
        
        task = loaddatastruct(animalinfo{2},animalinfo{3},'task');
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos');
        
        % retrieve epochs here
        epochs = evaluatefilter(task,epochmaker(epochstring));
         numeps = size(epochs,1);
         
        clear A
        A.animalname = animalname;
        A.epochs = epochs;
        A.corrwin = corrwin;
        A.kernsigs = kernsigs;
        A.epochstring = epochstring;
        A.statespec = statespec;
        A.meancorr = nan(numeps,numkernsigs,numstates);   %  [ epochnum x kernsignum x state]
        A.percent_nonsig = nan(numeps,numkernsigs,numstates); 
        
        for ee = 1:numeps
            
            d = epochs(ee,1);
            ep = epochs(ee,2);
            
            % load power traces
            trace = cell(1,2);
            timevec = cell(1,2);
            Fs = [nan nan];
            for g = 1:2
                out = loadtracestruct(animalinfo{2}, animalinfo{3},tracespec{g}{1}, d, ep);
                TF = tracespec{g}{2};
                trace{g} = zscorer(out{d}{ep}{TF}.powertrace);
                timevec{g} = out{d}{ep}{TF}.eegtimesvec_ref;
                Fs(g) = round(out{d}{ep}{TF}.samprate);
            end
            if ~all(timevec{1} == timevec{2})
                keyboard
            end
            
            for kk = 1:length(kernsigs)
               
                s = kernsigs(kk);
                sig = Fs(1) * s / 1000;
                kern = gaussian(sig,8*sig);
                trace1 = smoothvect(trace{1},kern);
                trace2 = smoothvect(trace{2},kern);
                
                numwins = floor(  (timevec{1}(end) - timevec{1}(1)) / corrwin ) ;
                corrtrace = nan(1,numwins);
                timetrace = nan(1,numwins);
                
                for w = 1:numwins
                    j = floor(   (w-1) * corrwin * Fs(1) + 1   );
                    k = floor(    w * corrwin * Fs(1)              );
                    m = round(mean([j k]));
                    timetrace(w) = timevec{1}(m);
                    [RHO,PVAL] = corr(trace1(j:k)',trace2(j:k)');
                    if PVAL < 0.05
                        corrtrace(w) = RHO;
                    end
                end
                
                % calculate correlation values in each state
                for state = 1:length(statespec)
                     [out, ~] = evaluatetimefilter(animalinfo{2},animalinfo{3},statespec{state},[d ep]);
                     includeperiods = out{d}{ep};
                     
                     goodinds = logical(isExcluded(timetrace,includeperiods));
                     
                     vals = corrtrace(goodinds);
                     meancorr = nanmean(vals);
                     
                     A.meancorr(ee,kk,state) = meancorr;   %  [ epochnum x kernsignum x state]
                    
                     % report percentage of windows with non-sig corr
                     percent_nonsig = 100 * sum(isnan(vals))/length(vals);
                        disp(num2str(percent_nonsig))
                        
                     A.percent_nonsig(ee,kk,state) = percent_nonsig;
                     
                     if plot_test
                         % sanity check
                         figure;
                         tvec = timevec{1} - timevec{1}(1);
                         plot(tvec, trace1 ,'-b','linewidth',2); hold on
                         plot(tvec, trace2 ,'-r','linewidth',2); hold on
                         yshift = -2;
                         tvec2 = timetrace - timevec{1}(1);
                         plot(tvec2, corrtrace + yshift,'ok','linewidth',2); hold on
                         plot(tvec2, corrtrace + yshift,'-k','linewidth',2); hold on
                         plot([tvec2(1) tvec2(end)]  , [-1 -1] + yshift,'k--','linewidth',1)
                         plot([tvec2(1) tvec2(end)]  ,  [1 1] + yshift,'k--','linewidth',1)
                         xlim([tvec(1) tvec(end)])
                         keyboard
                     end
                     
                end
                
            end
            
        end
        
        cd(datadir)
        save(sprintf('Gamma_runningcorr_%s',animalname(1:3)),'A')
        
    end
    

end



if plot_animal
    
    
    
    
    numanimals = length(animals_toplot);
    
    if numanimals > 1
        figure
        %H = figure('units','normalized','outerposition',[0 0 1 1]);
        fsize = 12;
    else
        H = figure('units','normalized','outerposition',[0 0 .8 .3]); 
        fsize = 16;
    end
    
    plotcounter = 1;
    
    for aa = 1:numanimals
        
        animalname = animals_toplot{aa};
        
        cd(datadir)
        load(sprintf('Gamma_runningcorr_%s.mat',animalname(1:3)),'A')
            numstates = length(A.statespec);
        
        for st = 1:numstates
           
            subplot(numanimals,numstates,plotcounter)
            
            if plot_allepochs
                for ep = 1:size(A.meancorr,1)
                    clr = [rand rand rand];
                    plot(A.kernsigs,A.meancorr(ep,:,st),'-','Color',clr,'linewidth',2); hold on
                    %pause
                    %close all
                end
            end
            
            naneps = isnan(prod(A.meancorr(:,:,st),2));
            meanplot = mean(A.meancorr(~naneps,:,st),1);
            semplot = std(A.meancorr(~naneps,:,st),0,1) / size(A.meancorr(~naneps,:,st),1);
            
            jbfill(A.kernsigs,meanplot+semplot,meanplot-semplot,[0 0 0],[0 0 0],0,1);
            
            title(sprintf('%s, %d epochs',animalname(1:3),sum(~naneps)),'fontsize',12,'fontweight','bold')
            if numanimals == 1
                xlabel('Kernel sigma (ms)','fontsize',fsize)
                ylabel('Average correlation (r)','fontsize',fsize)
            end
            ylim([-.2 .5])
            if plot_xlimit(2) == 500
                ylim([-.1 .2])
                set(gca,'ytick',[-.1 0 .2])
            else
                ylim([-.2 .5])
                set(gca,'ytick',[-.2 0 .5])
            end
            
            set(gca,'fontsize',fsize)            
            
            if ~isempty(plot_xlimit)
                xlim(plot_xlimit);
            end
            
            plotcounter = plotcounter + 1;
            
            
        end
        
    end
    
end




