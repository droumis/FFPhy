% Cluster decoding on W-track -- psuedo 1D "armdists" (Wu--Foster-2014 Fig 1)


datadir_pseudo = '/opt/data13/kkay/Superclustdecode_data';

calculate = 1;   % identifies candidate events and saves in P variable (Wu-Foster-2014)
plotFR = 1;    % forward vs. reverse comparison
plotLR = 0;    % local vs. remote

    if calculate
        targetfile = '';
       nn = 1;  % 1: local, 2: remote
       animals_tocalc = {'Bond'}; % leave blank if want to process all files
             animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       tracekern_ms = [];   
       clear tracespec
       tracespec{1} = {'rippletrace',1,[.7 .7 .7]} ;   % tracefile name + TFgroup
       tracespec{2} = {'fastgammatrace',2,[.4 .4 1]} ; 
       tracespec{3} = {'lowgammatrace',1,[1 0 0 ]}   ;
      
       tracespec{4} = {'fastgammatrace',1,[0 0 1]} ; 
       tracespec{5} = {'lowgammatrace',2,[1 .8 .2]} ; 
       
            numbands = length(tracespec);
    end
    
  

if calculate
    
    % initialize outputs
    DATA = cell(1,numbands);  % {tracenum}[  an   d   starttime   endtime   dur   power ]
    
    % load all files
    cd(datadir_pseudo)
    if isempty(targetfile)
        filenames = dir('*.mat');
    else
        filenames = dir(targetfile);
    end
        
    % load each file 
    for rr = 1:length(filenames)
       
        % load pseudo file  'P'
        cd(datadir_pseudo)
        filename = filenames(rr).name;
        load(filename,'P')

        animalname = P.animalname;
            animalinfo = animaldef(animalname);
            an = find(strcmp(animalname,animals_order));
            
        if ~any(strcmp(animalname,animals_tocalc))
            continue
        end
            
        d = P.dayep(1);
        ep = P.dayep(2);
        
        disp(sprintf('%s %d %d replay gamma',animalname(1:3),d,ep))

        tdel = P.binvec_c(2) - P.binvec_c(1);  % time bin size
        
        % load pos data for epstart
        pos = loaddatastruct(animalinfo{2}, animalinfo{3},'pos', d);
            epstart = pos{d}{ep}.data(1,1);
       
        outputblock = [];
        
        if ~isempty(P.replay{nn})

              starttimes = P.replay{nn}(:,1) - tdel/2;
              endtimes = P.replay{nn}(:,2) + tdel/2;
              frflags = P.replay{nn}(:,end);
                    numreplays = length(starttimes);
                    durations = endtimes - starttimes;
                    
              % get powers
              tracepowers = nan(numreplays,1);
              for g = 1:numbands
                  out = loadtracestruct(animalinfo{2}, animalinfo{3},tracespec{g}{1}, d, ep);
                  TF = tracespec{g}{2};
                    [ztrace,~,~] = zscoretrace(out{d}{ep}{TF},tracekern_ms);
                  timevec = out{d}{ep}{TF}.eegtimesvec_ref;
                  for rr = 1:numreplays
                      a = lookup(starttimes(rr),timevec);
                      b = lookup(endtimes(rr),timevec);
                      tracepowers(rr) = mean(ztrace(a:b));
                  end
                 outputblock = [ an * ones(numreplays,1) ...
                              d  * ones(numreplays,1) ...
                              ep * ones(numreplays,1) ...
                              starttimes ...
                              endtimes ...
                              durations ...
                              frflags ...
                              tracepowers ];                  
                  
                  DATA{g} = [ DATA{g}  ; outputblock ];
                  
              end
             

        end
        
    end
    
end



  if plotFR
        numbands = length(DATA);
        finds = DATA{1}(:,end-1) == 1 ;
        rinds = DATA{1}(:,end-1) == -1 ;
        
        fmean = []; fsem = []; numfor = sum(finds);
        rmean = []; rsem = []; numrev = sum(rinds);
        
        for bb = 1:numbands
            fvals = DATA{bb}(finds,end);
            rvals = DATA{bb}(rinds,end);
            fmean(bb) = mean(fvals);
            fsem(bb) = std(fvals)/sqrt(length(fvals));
            rmean(bb) = mean(rvals);
            rsem(bb) = std(rvals)/sqrt(length(rvals));
        end
        
        figure;
        for b = 1:numbands
            
            subplot(1,2,1)
            bar(b,fmean(b),'facecolor',tracespec{b}{3}); hold on
            errorbar(b,fmean(b),fsem(b),'.k','linewidth',1)
            title(sprintf('Forward replay (n = %d) (kern: %d ms)',numfor,tracekern_ms),'fontsize',11,'fontweight','bold')
            set(gca,'fontsize',12)
            ylabel('Mean Z-score power')
            set(gca,'xtick',[])
            xlabel('Gamma type')
            ylim([-1 2]);
            
            subplot(1,2,2)
            bar(b,rmean(b),'facecolor',tracespec{b}{3}); hold on
            errorbar(b,rmean(b),rsem(b),'.k','linewidth',1)
            title(sprintf('Reverse replay (n = %d)',numrev),'fontsize',11,'fontweight','bold')
            xlabel('Gamma type')
             ylabel('Mean Z-score power')
             set(gca,'xtick',[])
            set(gca,'fontsize',12)
            ylim([-1 2]);
        end
  end

    
    if plotLR
        numbands = length(DATA);
        finds = DATA{1}(:,end-1) == 1 ;
        rinds = DATA{1}(:,end-1) == -1 ;
        
        fmean = []; fsem = []; numfor = sum(finds);
        rmean = []; rsem = []; numrev = sum(rinds);
        
        for bb = 1:numbands
            fvals = DATA{bb}(finds,end);
            rvals = DATA{bb}(rinds,end);
            fmean(bb) = mean(fvals);
            fsem(bb) = std(fvals)/sqrt(length(fvals));
            rmean(bb) = mean(rvals);
            rsem(bb) = std(rvals)/sqrt(length(rvals));
        end
        
        figure;
        for dd = 1:numbands

            %forward
            subplot(1,2,1)
            bar(dd,fmean(dd)); hold on
            errorbar(dd,fmean(dd),fsem(dd),'linewidth','none')
            title(sprintf('forward replay (n = %d) (kern: %d ms)',numfor,tracekern_ms))
            ylim([-1 2]);

            % reverse
            subplot(1,2,2)
            bar(dd,rmean(dd)); hold on
            errorbar(dd,rmean(dd),rsem(dd),'linespec','k','linewidth','none')
            title(sprintf('reverse replay (n = %d)',numrev))
            ylim([-1 2]);
            
        end
        
    end
