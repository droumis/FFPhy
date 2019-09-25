function sss_sdss2arcplots(nfwaves,nftimes,cls,unqcs,nsweeps,sweepd,stimonset,window,arfp,plot_ampl,plot_sd, plot_type)

% SDSS2ARCPLOTS is a plotting function called by sdss2arc. See SDSS2ARC for the purpose of the program. 
% All rights reserved
% Tansu Celikel
% plot the raw waveforms with the firing statistics

% MODIFIED: SHANTANU JADHAV (FEB 2005)


%stimonset=500; window=499;

ctr=1;
for lpfs=1:ceil(length(unqcs)/5)
    figure ([99+lpfs])
    
    if length(unqcs) >= ctr+4 
        tempcls=unqcs(ctr:ctr+4);
    else
        tempcls=unqcs(ctr:end);
    end
    
    tctr=1;
    for lptempcls=1:length(tempcls)
        
        sp_idx=find (cls==tempcls(lptempcls));
        
         m = mean(nfwaves(sp_idx,:));
         s = std(nfwaves(sp_idx,:));
        % raw waveform
        subplot (5,5,tctr); hold on;   %subplot (length(tempcls),5,tctr);
        if plot_type==0
            plot (nfwaves(sp_idx,:)', 'k'); axis tight
        else
            [n,x,y] = hist2d(nfwaves);
            imagesc(x,y,n); axis xy; colormap jet;
        end
        plot(m, 'c', 'LineWidth', 2);
        if plot_sd==1
            plot(m+2*s, 'r.', 'LineWidth', 2);
            plot(m-2*s, 'r.', 'LineWidth', 2);
        end
        axis([0 size(nfwaves,2) -plot_ampl plot_ampl])
        if lptempcls ==1; title ('Raw Waves');end 
        ylabel (['CL- ' num2str(tempcls(lptempcls))])
        dtimes=diff(nftimes(sp_idx)/1000);
        yaxisval=get(gca,'Ylim');
        text (1,yaxisval(1,2)*0.6,['Sp<ISI=' num2str(length(find ((dtimes >= 0) & (dtimes <= arfp/1000))))]);
        text (1,yaxisval(1,2)*0.8,['Nsp=' num2str(length(sp_idx))]);
        
        % raster plot
        subplot (5,5,tctr+1);
        sptms=mod(nftimes(sp_idx),10000);
        sptrls=(nftimes(sp_idx)-sptms)/10000;
        plot (round(sptms/10),sptrls, 'b.');
        set(gca, 'Xlim', [0 sweepd])
        
        % PSTH -full sweep
        subplot (5,5,tctr+2);    
        [xvals]=histc(round(sptms/10),1:1:sweepd);
        stairs (xvals); axis tight
        
        % PSTH - stim onset with background blown up
        subplot (5,5,tctr+3);
        stairs (xvals(stimonset-window:stimonset+window)); axis tight;
        if lptempcls ==1; title (['Stimonset +/- ' num2str(window) 'ms']);end 
        
        % Spike count
        subplot (5,5,tctr+4);
        bcksp=sum(xvals(stimonset-window:stimonset))/(nsweeps);
        evosp=sum(xvals(stimonset:stimonset+window))/(nsweeps);
        bar ([1;2],[bcksp;evosp]); axis tight;
        if lptempcls ==1; title ('Spikes per stimulus');end 
        if lptempcls==length(tempcls);      
            set (gca, 'XTickLabel', {'Background' 'StimR'})
        else
            set (gca, 'XTickLabel', [])
        end
        tctr=tctr+5;
        
    end
    ctr=ctr+5; clear tctr
    
   redimscreen100(99+lpfs)
    
end
