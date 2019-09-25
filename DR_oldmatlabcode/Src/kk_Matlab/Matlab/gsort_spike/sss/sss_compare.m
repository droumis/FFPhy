function sss_compare (spikes,clustnum,arfp,Coloropt, parameters)

waveforms = spikes.waveforms;
times = spikes.fstimes; 
newplots = spikes.hierarchy.assigns;

store_times=[];

    clr = ['w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r' 'w' 'y' 'm' 'b' 'c' 'r'];
        
    h=figure; redimscreen_hor; whitebg(h)
    Tsp=zeros(size(times,1),1);
    
    for lp=1:size(clustnum,1)
        sp_ind=(find(newplots==clustnum(lp)));
        
        if (Coloropt==0); % if units combined 
            temptimes=times(sp_ind);
            tempwaves=waveforms(sp_ind,:);          
            store_times=[store_times;temptimes];
            
            subplot (4,6,[1:18])
            plot(tempwaves',clr(1)); hold on;
            plot (mean (tempwaves), clr(4), 'Linewidth', 2);
            %set (gca, 'XLim',[1 size(spikes.waveforms,2)]); 
            %yl=get(gca, 'Ylim');
            
            if lp==size(clustnum,1)
                set (gca, 'XLim',[1 size(spikes.waveforms,2)]); 
                yl=get(gca, 'Ylim');
                dtimes=diff(sort(store_times/1000));        
                title(['Clusters ' num2str(clustnum') ' | Nspikes=  ' num2str(length(store_times)) ' | <ISI spikes=  ' num2str(length(find ((dtimes<arfp/1000) & (dtimes >0))))])
            end
            
            subplot (4,6,[19:24])
            axis off;
            text(0,1+(-0.2*lp),['Cluster ' num2str(clustnum(lp)) ' (' clr(lp) ')']);
            text(0.2,1+(-0.2*lp),['Nspikes=  ' num2str(length(temptimes))]);
            dtimes=diff(sort(temptimes/1000));  
            text(0.5,1+(-0.2*lp),['<ISI spikes=  ' num2str(length(find ((dtimes<arfp/1000) & (dtimes >0))))]); 
            
            
            
        else  % if units are NOT combined
            temptimes=times(sp_ind);
            tempwaves=waveforms(sp_ind,:);
            store_times=[store_times;temptimes];
            
            subplot (4,6,[1:18]) 
            plot(tempwaves',clr(lp)); hold on;   
            plot (mean (tempwaves), clr(lp+3), 'Linewidth', 2);
            %set (gca, 'XLim',[1 size(spikes.waveforms,2)]); 
            %ylabel ('Voltage (uV)');
            
            if lp==size(clustnum,1)
                set (gca, 'XLim',[1 size(spikes.waveforms,2)]); 
                yl=get(gca, 'Ylim');
                dtimes=diff(sort(store_times/1000));        
                title(['Clusters ' num2str(clustnum') ' | Nspikes=  ' num2str(length(store_times)) ' | <ISI spikes=  ' num2str(length(find ((dtimes<arfp/1000) & (dtimes >0))))])
            end
            
            subplot (4,6,[19:24])
            axis off;
            text(0,1+(-0.2*lp),['Cluster ' num2str(clustnum(lp)) ' (' clr(lp) ')']);
            text(0.2,1+(-0.2*lp),['Nspikes=  ' num2str(length(temptimes))]);
            dtimes=diff(sort(temptimes/1000));  
            text(0.5,1+(-0.2*lp),['<ISI spikes=  ' num2str(length(find ((dtimes<arfp/1000) & (dtimes >0))))]); 
            
           
        end
            
        if lp==1, t1=temptimes/1000; clust1=clustnum(lp); end
        if lp==2, t2=temptimes/1000; clust2=clustnum(lp); end
        
        clear sp_ind temptimes tempwaves;
    end
    
    energy=spikes.hierarchy.interface_energy(clust1,clust2);
    if (isfield(spikes,'parameters'))
        tmin = spikes.parameters.tmin; tref = spikes.parameters.tref;
    elseif (length(parameters.tmin)~=0)
        tmin = parameters.tmin; tref = parameters.tref;
    else
        tmin=0.001; tref=0.002;
    end

    [allow, scores] = isiQuality(t1, t2, tmin, 0.010, tref, spikes.Fs);
    subplot (4,6,[19:24])
    text(0,0.5+(-0.2*(lp)),['ISI Quality Allow = ' num2str(allow)]);
    text(0,0.2+(-0.2*(lp)),['ISI Scores: Clus 1st,2nd,Comb = ' num2str(scores(1)) '; ' num2str(scores(2)) '; ' num2str(scores(3))]);
    text(0.2,0.5+(-0.2*(lp)),['Energy = ' num2str(energy)]);

    

