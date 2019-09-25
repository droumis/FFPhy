function [spikes,aspkfile]=sss_joinsortedfiles_cont3(files)

%% sss_joinsortedfiles_cont3({'spike_cut_tet_thr5_allpart1a_sort1.mat'; 'spike_cut_tet_thr5_allpart1b_sort1rename.mat'; 'spike_cut_tet_thr5_allpart2_sort1rename.mat'});
%% eg. spec1='spike_cut_'; spec2='_tet_thr5',ch=4; nfiles=10;
%% pushfiles=0 for no, 1 for yes

%%%%% DO NOT STORE INDIVIDUAL CHANNEL WAVEFORMS: SAVE SPACE %%%%%%%

%%%% Initialize finalallspks

ch=1:4;
finalallspks.waveforms=[];

for nch=ch(1):ch(4)
    cmd=sprintf('finalallspks.waveforms_ch%d = [];',nch); eval(cmd);
end

finalallspks.spiketimes=[]; finalallspks.swtimes=[];
finalallspks.fstimes=[];finalallspks.ftimes=[];
finalallspks.sweep = [];finalallspks.trial = [];
finalallspks.stimulus = []; finalallspks.igorstim = [];finalallspks.nsweeps = 0;
finalallspks.hierarchy.assigns=[];
finalallspks.sweepall=[];

for n=1:size(files,1)
    n
    filename=[files{n}];
    load(filename);

    for nch=ch(1):ch(4),
        cmd=sprintf('finalallspks.waveforms_ch%d = [finalallspks.waveforms_ch%d;  spikes.waveforms_ch%d];',nch, nch, nch); eval(cmd);
    end
    
    
    finalallspks.spiketimes=[finalallspks.spiketimes;spikes.spiketimes]; finalallspks.ftimes=[finalallspks.ftimes;spikes.ftimes];
    finalallspks.swtimes=[finalallspks.swtimes;spikes.swtimes]; finalallspks.fstimes=[finalallspks.fstimes;spikes.fstimes];  %%% UPDATED FSTIMES
    finalallspks.sweep = [finalallspks.sweep; spikes.sweep]; finalallspks.trial = [finalallspks.trial; spikes.trial];
    finalallspks.stimulus = [finalallspks.stimulus; spikes.stimulus];
    finalallspks.igorstim = [finalallspks.igorstim; spikes.igorstim];
    finalallspks.nsweeps = size(files,1);
    
    finalallspks.sweepall=[finalallspks.sweepall, spikes.sweepall];
    
    finalallspks.hierarchy.assigns=[finalallspks.hierarchy.assigns; spikes.hierarchy.assigns];
    
    finalallspks.Fs=spikes.Fs; finalallspks.threshT=spikes.threshT;
    finalallspks.threshV=spikes.threshV;
    finalallspks.stimonset=spikes.stimonset; finalallspks.window=spikes.window; finalallspks.sweepd=spikes.sweepd;
   
    %if n==size(files,1), save temp finalallspks; end
    
    %finalallspks.waveforms=[finalallspks.waveforms; spikes.waveforms];
    clear spikes
    
end % file loop

clear n files 
clear ch
finalallspks.nsweeps=floor(max(finalallspks.fstimes)/1000);


%%%%% DO NOT COMBINE FOR FINAL WAVEFORMS: TOO BIG%%%%%%%%%%%%%%
% for nch=1:ch
%     cmd=sprintf('finalallspks.waveforms = [ finalallspks.waveforms; finalallspks.waveforms_ch%d];',nch); eval(cmd);
% end
%finalallspks.waveforms = finalallspks.waveforms;

%save ('finalfile', 'finalallspks');
spikes=finalallspks;
save ('finalfile2', 'spikes');


