function spikes = finalize_sort(spikes)

%--------------------------------------------------------------------------

warning off; 


tmin=0.001;         % TMIN WILL BE IN SPIKES.PARAMETERS IN FUTURE VERSIONS
arfp = tmin*1000;
clusters = spikes.hierarchy.assigns;
unqcs =  unique(clusters);

sss_sdss2arcplots(spikes.waveforms, spikes.fstimes ,clusters, unqcs, spikes.nsweeps, spikes.sweepd, spikes.stimonset, spikes.window, arfp);
pause
sss_compwaves2 (spikes.waveforms, spikes.fstimes, clusters, arfp);    %ARPF FOR ISI = 1MS;

%combine unit
combunits=input('If you want to combine clusters, please enter cluster numbers now (i.e [1 2 3])...');

while length(combunits)>=1
    
    for lpnc=2:length(combunits)
        spikes = merge_clusters(spikes, combunits(1), combunits(lpnc));
        clusters = spikes.hierarchy.assigns;
        unqcs =  unique(clusters);
    end
    
    close all
    sss_sdss2arcplots(spikes.waveforms, spikes.fstimes ,clusters, unqcs, spikes.nsweeps, spikes.sweepd, spikes.stimonset, spikes.window, arfp);
    pause
    sss_compwaves2 (spikes.waveforms, spikes.fstimes, spikes.hierarchy.assigns, arfp);    %ARPF FOR ISI = 1MS;
    
    % sdss2arcplots(nfwaves,nftimes,cls,unqcs,nsweeps,sweepd,stimonset,window,fig_opt,pt_opt,reminder)
    % compwaves2 (nfwaves, nftimes, newassigns,1);
    
    combunits= input('To CONTINUE reclustering, you can enter cluster numbers now (i.e [1 2 3])...');
end


% CLUSTERS FINALIZED, AND TREE UPDATED: NOW SORT ASSIGNMENTS 
% SORT ASSIGNMENTS:  SEND TO SPIKES.NEWASSIGNS, AND KEEP MAP
% SO SPIKES.HIERARCHY.ASSIGNS ALWAYS HAS A FORM COINCIDING WITH THE TREE

% ONLY THING NOT UPDATED IS SPIKES.CLUSTERS 

[spikes.newassigns, spikes.clustermap] = sortAssignments2(spikes.hierarchy.assigns);
clusters = spikes.newassigns;
unqcs =  unique(clusters);


save tempsdss