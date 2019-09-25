% use this script to visually check theta epochs extracted by
% gethighthetatimes -- also use to empirically determine good values for
% powerratio

% 1st block: checks for BIMODALITY
% 2nd block: plots eeg, theta, powerratios, velocity for a chosen tetrode
% 3rd block: plots high theta times for all tetrodes in animal + ripple
        % times


animaldir = '/data12/mari/Egy/';
animalprefix = 'egy';
epochs = [6 8];
tetlist = [];

% plot options
regions = {'cc' 'ctx' 'Reference' 'CA1' 'CA2' 'CA3'};
rem_flag = 1;
nrem_flag = 1;
spindle_flag = 1;
scroll_plot = 0;
event_plot = 'rem';
    event_tet = 15;
  % mua parameters
kernelsize = .010;          % in s
muaregions = {'ctx'};   % regions for which to plot mua
multi_ctxconsensus_flag = 1;  % plots ctx-consensus mua
multi_smooth_flag = 0;        % plots smoothed mua for same tetrode
multi_raster_flag = 0;        % plot raster of mua for same tetrode
    

% set these to taste
velocitythresh = 5;
windowsize = 10;            % in seconds
mindur = 0.5;
ntheta_thresh = 3;
consensus = 0;            % set to 1 (see below) if want to plot same acceptable period for all
                          % set to 0 if plot each individual tetrode's
                          % acceptable period


pos = loaddatastruct(animaldir, animalprefix, 'pos', epochs(1));
    pos = pos{epochs(1)}{epochs(2)};
    windowsize_possamp = windowsize*29.97;
    windowsize_eegsamp = windowsize*1500;
    nowindows = floor((pos.data(end,1)-pos.data(1,1))/windowsize);
    Fs = 1500;
    timevec = 0:(1/Fs):windowsize;
tetlist = [1:30];
tetinfo = loaddatastruct(animaldir, animalprefix,'tetinfo');

ripples = loaddatastruct(animaldir, animalprefix,'ripples',epochs(1));
sleep = loaddatastruct(animaldir, animalprefix,'sleep',epochs(1));
rem = loaddatastruct(animaldir, animalprefix,'rem',epochs(1));
nrem = loaddatastruct(animaldir, animalprefix,'nrem',epochs(1));
multi = loaddatastruct(animaldir, animalprefix,'multi',epochs(1));
spindles = loaddatastruct(animaldir, animalprefix,'spindles',epochs(1));

    
    % retrieve recording areas of each tetrode and determine whether valid to plot
    validtet = zeros(size(tetlist));
    hpctet = zeros(size(tetlist));
    ctxtet = zeros(size(tetlist));
    
    tetlist_region = cell(size(tetlist));
    for tet=1:length(tetlist)
        try
            tetlist_region{tet} = tetinfo{epochs(1)}{epochs(2)}{tetlist(tet)}.area;
            % plot Reference, cc, and ctx no matter what
            if strcmp(tetlist_region{tet},'Reference') || strcmp(tetlist_region{tet},'cc') || strcmp(tetlist_region{tet},'ctx')
                validtet(tet) = 1;
                disp(sprintf('tetrode %d : reference',tetlist(tet)))
                if strcmp(tetlist_region{tet},'ctx')
                    ctxtet(tet) = 1;
                end
                continue
            end
            % flag hpc tetrodes
            if strcmp(tetlist_region{tet},'CA1') || strcmp(tetlist_region{tet},'CA2') || strcmp(tetlist_region{tet},'CA3') || strcmp(tetlist_region{tet},'DG')
                hpctet(tet) = 1;
            end
            % plot other tetrodes only if two clustered cells
            try 
                nocellsontet{tetlist(tet)} = tetinfo{epochs(1)}{epochs(2)}{tetlist(tet)}.numcells;
                if nocellsontet{tetlist(tet)} >= 2
                    validtet(tet) = 1;
                    disp(sprintf('tetrode %d : %d cells',tetlist(tet),nocellsontet{tetlist(tet)}))
                    continue
                else
                    disp(sprintf('tetrode %d : ignored, %d cells',tetlist(tet),nocellsontet{tetlist(tet)}))
                    continue
                end
            catch
                disp(sprintf('tetrode %d : no .numcells field',tetlist(tet)))
                continue
            end
        catch
            disp(sprintf('tetrode %d : no area field / other error',tetlist(tet)))
            continue
        end
    end
    
    
    


    
    
    
    
    dumflag = 1;
    
    % load full eeg, sleep, rem, ripple, nrem, spindles times
    for tet=1:length(tetlist)       

        if validtet(tet)
            
            % eeg
            eeg{tet} = loadeegstruct(animaldir, animalprefix, 'eeg', epochs(1), epochs(2),tetlist(tet));
            eeg{tet} = eeg{tet}{epochs(1)}{epochs(2)}{tetlist(tet)};
            eegtimes = geteegtimes(eeg{tet})';
            

            
            % sleep times (only do once)
            if dumflag
                sleeptrace = nan(size(eegtimes))';
                if ~isempty(sleep{epochs(1)}{epochs(2)})
                    sleepduration = round(sum(sleep{epochs(1)}{epochs(2)}.endtime - sleep{epochs(1)}{epochs(2)}.starttime));
                    numsleep = length(sleep{epochs(1)}{epochs(2)}.starttime);
                    for j=1:numsleep
                        startind = lookup(sleep{epochs(1)}{epochs(2)}.starttime(j),eegtimes);
                        endind = lookup(sleep{epochs(1)}{epochs(2)}.endtime(j),eegtimes);
                        sleeptrace(startind:endind) = 1;
                    end
                else
                    sleepduration = 0;
                end
                disp(sprintf('Total sleep this epoch: %d s',round(sleepduration)))
                dumflag = 0;
            end
            
            % ripple times
            if hpctet(tet)
            rippletrace{tet} = nan(size(eegtimes))';
            noripples = length(ripples{epochs(1)}{epochs(2)}{tetlist(tet)}.startind);
            %disp(sprintf('tetrode %d : %d cells, %d ripples',tetlist(tet),nocellsontet{tetlist(tet)},noripples))
            for j=1:noripples
                startind = ripples{epochs(1)}{epochs(2)}{tetlist(tet)}.startind(j);
                endind = ripples{epochs(1)}{epochs(2)}{tetlist(tet)}.endind(j);
                rippletrace{tet}(startind:endind) = 1;
            end
            end
            
            % rem times
            if hpctet(tet)
            remtrace{tet} = nan(size(eegtimes))';
            if ~isempty(rem{epochs(1)}{epochs(2)}{tetlist(tet)})
                remduration = round(sum(rem{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime - rem{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime));
                norem = length(rem{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime);
                for j=1:norem
                    startind = lookup(rem{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime(j),eegtimes);
                    endind = lookup(rem{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime(j),eegtimes);
                    remtrace{tet}(startind:endind) = 1;
                end
            else
                remduration = 0;
            end
            end
            
            % nrem times
            nremtrace{tet} = nan(size(eegtimes))';
            if ~isempty(nrem{epochs(1)}{epochs(2)}{tetlist(tet)})
                nremduration = round(sum(nrem{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime - nrem{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime));
                nonrem = length(nrem{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime);
                for j=1:nonrem
                    startind = lookup(nrem{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime(j),eegtimes);
                    endind = lookup(nrem{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime(j),eegtimes);
                    nremtrace{tet}(startind:endind) = 1;
                end
            else
                nremduration = 0;
            end
            
            % spindle times
            spindletrace{tet} = nan(size(eegtimes))';
            if ~isempty(spindles{epochs(1)}{epochs(2)}{tetlist(tet)})
                spindleduration = round(sum(spindles{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime - spindles{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime));
                numspindles = length(spindles{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime);
                for j=1:numspindles
                    startind = lookup(spindles{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime(j),eegtimes);
                    endind = lookup(spindles{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime(j),eegtimes);
                    spindletrace{tet}(startind:endind) = 1;
                end
            else
                spindleduration = 0;
            end  
            
            % report durations
            try
                disp(sprintf('tetrode %d : %d cells, %d s of REM',tetlist(tet),nocellsontet{tetlist(tet)},remduration))
            catch
                disp(sprintf('tetrode %d : 0 cells, %d s of REM',tetlist(tet),remduration))
            end
            try
                disp(sprintf('tetrode %d : %d cells, %d s of NREM',tetlist(tet),nocellsontet{tetlist(tet)},nremduration))
            catch
                disp(sprintf('tetrode %d : 0 cells, %d s of NREM',tetlist(tet),nremduration))
            end            
            try
                disp(sprintf('tetrode %d : %d cells, %d s of Spindles',tetlist(tet),nocellsontet{tetlist(tet)},spindleduration))
            catch
                disp(sprintf('tetrode %d : 0 cells, %d s of Spindles',tetlist(tet),spindleduration))
            end                       
        end
    end

    
    % collect all MUA activity from ctx tetrodes
    ctxmua = [];
    dummyflag = 0;
    kernel = gaussian(kernelsize*Fs,8*kernelsize*Fs);
            
    %figure
    %hold on
    for tet = 1:length(tetlist)        
        if ~isempty(tetlist_region{tet})
            if ~dummyflag
                ctxmuatimes = geteegtimes(eeg{tet})';
                ctxreftet = tet;
                dummyflag = 1;
            end
            mua = multi{epochs(1)}{epochs(2)}{tet}/10000;
            if strcmp(tetlist_region{tet},'ctx')
                ctxmua = [ctxmua ; mua];
            end
            % bin the spikes
            N = hist(mua,ctxmuatimes);
            % smooth the rate
            firing{tet} = smoothvect(N,kernel) * Fs;
            %plot(eegtimes-eegtimes(1),firing,'Color',[rand rand rand])
        end
    end
    
    % combine mua from all all cortical tetrodes
    N_all = hist(ctxmua,ctxmuatimes);
    firing_all = ( smoothvect(N_all,kernel) * Fs ) / sum(ctxtet);
    
    
    
    
    if event_plot
        
        eval(sprintf('events = %s ;',event_plot));
        % plot events from manually chosen event_tet OR
            % tet w/ longest duration of events :
        if isempty(event_tet)
            event_tet = [];
            dummy = 0;
            for tet = tetlist
                totaleventdur = sum(events{epochs(1)}{epochs(2)}{tetlist(tet)}.endtime - events{epochs(1)}{epochs(2)}{tetlist(tet)}.starttime);
                if totaleventdur > dummy
                    event_tet = tet;
                end
            end
        end
        
        eventdata = events{epochs(1)}{epochs(2)}{tetlist(event_tet)};
        numevents = length(eventdata.starttime);
        
        for w = 1:numevents
            
            midtime = (eventdata.endtime(w) + eventdata.starttime(w))/2;
            starttime = midtime - windowsize/2 ;
            endtime = midtime + windowsize/2 ;
            
            h = figure;
            %figure('units','normalized','outerposition',[0 0 1 1])
            
            %pos
            subplot(6,1,6)
            posindices = (lookup(starttime,pos.data(:,1))):(lookup(endtime,pos.data(:,1)));
            velocity = pos.data(posindices,5);
            threshindices = velocity > velocitythresh;
            hold on
            plot(velocity,'Color',[.8 .8 .8],'LineWidth',5)
            plot(velocity.*threshindices,'k.','MarkerSize',15)
            plot(velocitythresh*ones(1,length(velocity)),'k--')
            title('velocity')
            ylim([0 20])
            axis tight;
            
            %raweeg
            subplot(6,1,1:5)
            hold on
            
            plotcount=0;

            tetrodes_plotted = [];
            
            for r=1:6     % iterate through regions
                region=regions{r};
                for tet=[find(validtet==1)]
                    if strcmp(tetlist_region{tet},region)
                        
                        tetrodes_plotted = [tetrodes_plotted tet];
                        
                        eegtimes = geteegtimes(eeg{tet});
                        eegstartindex = lookup(starttime,eegtimes);
                        ctxmua_startindex=lookup(starttime,ctxmuatimes);
                        eegtrace = eeg{tet}.data(eegstartindex:(eegstartindex+windowsize_eegsamp))';
                        
                        %plot
                        %select colors
                        switch region
                            case 'CA1'
                                clr1 = [.9 .9 .9];
                                clr2 = [0 0 0];
                            case 'CA3'
                                clr1 = [1 .9 .9];
                                clr2 = [1 0 0];
                            case 'CA2'
                                clr2 = [0 .7 0];
                                clr1 = [.7 1 .7];
                            case 'Reference'
                                clr1 = [1 .9 1];
                                clr2 = [1 0 1];
                            case 'cc'
                                clr1 = [1 .9 1];
                                clr2 = [1 0 1];
                            case 'ctx'
                                clr1 = [1 .9 1];
                                clr2 = [1 0 1];
                        end
                     

                        
                        %eeg
                        plot(timevec,eegtrace-plotcount*1000,'Color',clr1,'LineWidth',1.2)
                        %sleep periods
                        plot(timevec,eegtrace.*sleeptrace(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000,'Color',clr2,'LineWidth',1.2)
                        
                        %spindletrace
                        if spindle_flag
                            if strcmp(region,'ctx')
                                 % bar format
                            plot(timevec,spindletrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'Color',[102 153 204]/255,'LineWidth',4)
                                % thick trace format   
                            plot(timevec,eegtrace.*spindletrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000,'Color',[102 153 204]/255,'LineWidth',1)
                            end
                        end

                        % multi
                        if multi_ctxconsensus_flag || multi_smooth_flag
                            if any(strcmp(region,muaregions))
                                if multi_ctxconsensus_flag
                                    muatrace = firing_all(ctxmua_startindex:(ctxmua_startindex+windowsize_eegsamp));
                                else
                                    muatrace = firing{tet}(ctxmua_startindex:(ctxmua_startindex+windowsize_eegsamp));
                                end
                                plotratio = range(eegtrace) / range(muatrace) ;
                                muatrace = plotratio * muatrace + min(eegtrace-plotcount*1000);
                                plot(timevec,muatrace,'Color',[.8 .8 .8],'LineWidth',1)
                            end
                        end
                        if multi_raster_flag
                            if any(strcmp(region,muaregions))
                                muaspikes = multi{epochs(1)}{epochs(2)}{tet}/10000;
                                    muaspikes = (muaspikes(muaspikes > starttime & muaspikes < endtime) - starttime);
                                    for s = 1:length(muaspikes)
                                        p = patchline([muaspikes(s) muaspikes(s)],[-plotcount*1000-300 -plotcount*1000+300],'linestyle','-','edgecolor','k','linewidth',1,'edgealpha',0.2);
                                    end
                                     %p(1) = patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
                                    %p(2) = patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
                                %line([muaspikes muaspikes],[-plotcount*1000-300 -plotcount*1000+300],'alpha',0.5,'Color',[.85 .85 .85],'LineWidth',1)
                            end          
                        end
                        
                        %rippletrace
                        if hpctet(tet)
                        plot(timevec,rippletrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'k','LineWidth',7)
                        end
                        
                        %remtrace
                        if hpctet(tet)
                        if rem_flag
                                % bar format
                            plot(timevec,remtrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'Color',[155 48 255]/255,'LineWidth',4)
                                % thick trace format   
                            plot(timevec,eegtrace.*remtrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000,'Color',[155 48 255]/255,'LineWidth',1)
                        end
                        end
                        
                        %nremtrace
                        if nrem_flag
                            plot(timevec,nremtrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'Color',[.7 .7 .7],'LineWidth',2)
                        end

                        plotcount=plotcount+1;
                    
                    end
                end
            end
            ylim([-plotcount*1000 400])
            line1 = sprintf('%s sleep day %d, epoch %d',animalprefix,epochs(1),epochs(2));
            line3 = sprintf('tetrodes: %s',mat2str(tetrodes_plotted));
            title({line1 line3},'FontSize',14,'FontWeight','bold')
            pause
            close all
        end
        
    end
        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        if scroll_plot
        for w = 1:nowindows
            starttime = pos.data(1,1) + (w-1)*windowsize;
            endtime = starttime + windowsize;
            h = figure;
            figure('units','normalized','outerposition',[0 0 1 1])
            
            %pos
            subplot(6,1,6)
            posindices = (lookup(starttime,pos.data(:,1))):(lookup(endtime,pos.data(:,1)));
            velocity = pos.data(posindices,5);
            threshindices = velocity > velocitythresh;
            hold on
            plot(velocity,'Color',[.8 .8 .8],'LineWidth',5)
            plot(velocity.*threshindices,'k.','MarkerSize',15)
            plot(velocitythresh*ones(1,length(velocity)),'k--')
            title('velocity')
            ylim([0 20])
            axis tight;
            
            %raweeg
            subplot(6,1,1:5)
            hold on
            
            plotcount=0;
            regions = {'cc' 'ctx' 'Reference' 'CA1' 'CA2' 'CA3'};
            tetrodes_plotted = [];
            
            for r=1:6     % iterate through regions
                region=regions{r};
                for tet=[find(validtet==1)]
                    if strcmp(tetlist_region{tet},region)
                        tetrodes_plotted = [tetrodes_plotted tet];
                        eegtimes = geteegtimes(eeg{tet});
                        eegstartindex = lookup(starttime,eegtimes);
                        eegtrace = eeg{tet}.data(eegstartindex:(eegstartindex+windowsize_eegsamp))';
                        %plot
                        %select colors
                        switch region
                            case 'CA1'
                                clr1 = [.9 .9 .9];
                                clr2 = [0 0 0];
                            case 'CA3'
                                clr1 = [1 .9 .9];
                                clr2 = [1 0 0];
                            case 'CA2'
                                clr1 = [0 .7 0];
                                clr2 = [.7 1 .7];
                            case 'Reference'
                                clr1 = [1 .9 1];
                                clr2 = [1 0 1];
                            case 'cc'
                                clr1 = [1 .9 1];
                                clr2 = [1 0 1];
                            case 'ctx'
                                clr1 = [1 .9 1];
                                clr2 = [1 0 1];
                        end
                        %eeg
                        plot(timevec,eegtrace-plotcount*1000,'Color',clr1,'LineWidth',1.2)
                        %highthetatrace
                        plot(timevec,highthetatrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000,'Color',clr2,'LineWidth',1.2)
                        %remtrace
                        if rem_flag
                            plot(timevec,remtrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'m','LineWidth',10)
                        end
                        if nrem_flag
                            plot(timevec,nremtrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'m','LineWidth',10)
                        end                        
                        %rippletrace
                        plot(timevec,rippletrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'b','LineWidth',5)
                        
                        plotcount=plotcount+1;
                    end
                end
            end
            ylim([-plotcount*1000 400])
            line1 = sprintf('%s highthetatimes day %d, epoch %d',animalprefix,epochs(1),epochs(2));
            line3 = sprintf('tetrodes: %s',mat2str(tetrodes_plotted));
            title({line1 line3},'FontSize',14,'FontWeight','bold')
            pause
            close all
        end
    end
    