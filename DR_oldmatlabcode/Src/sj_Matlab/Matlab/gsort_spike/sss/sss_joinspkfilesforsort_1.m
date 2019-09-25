function [spikes,aspkfile]=sss_joinspkfilesforsort_1(files)

%% sss_joinspkfilesforsort_1({'spiketetcut_all_thr-55_reduced8';'spiketetcut_all_thr-55_reduced8-2'});

%% Equivalent to sss_joinspikefile_cont5, but push times in a better way so
%% that it is easy to deconstruct spike times after sort to save in
%% individual files. Has to have Similar Spike Cut Parameters

%% Assume that no file will run longer than 10000secs/ safer/ use
%% 10^5+swtimes, and 10^8*fstimes for each file and push

%% For joining spk files in multiple directories before sorting. By
%% default, you have to adjust spike times to prevent overlap of spktimes
%% during sorting
%% POST Sort, deconstruct the spikesagain based on fstimes, and save in
%% indiv files: Use 


%%%%%%%%%%%%

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


push_swtimes=0; push_fstimes=0;

%%Start loop
for n=1:length(files)
    n
    filename=[files{n}];
    load(filename);
    
    if n==1,
        finalallspks.Fs=spikes.Fs;
        finalallspks.threshT=spikes.threshT;
        finalallspks.threshV=spikes.threshV;
        if isfield(spikes,'nsecs'), finalallspks.nsecs=0; end
        if isfield(spikes,'nsecs_reduced'); finalallspks.nsecs_reduced=0; end
    end

    finalallspks.spiketimes=[finalallspks.spiketimes; spikes.spiketimes]; 
    finalallspks.ftimes=[finalallspks.ftimes; spikes.ftimes ];               
    
    push_swtimes=n*10^5; push_fstimes=n*10^8;
    finalallspks.swtimes=[finalallspks.swtimes; spikes.swtimes + push_swtimes];              
    finalallspks.fstimes=[finalallspks.fstimes; spikes.fstimes + push_fstimes];      
    
    finalallspks.waveforms=[finalallspks.waveforms; spikes.waveforms];
    for nch=ch(1):ch(4),
        cmd=sprintf('finalallspks.waveforms_ch%d = [finalallspks.waveforms_ch%d;  spikes.waveforms_ch%d];',nch, nch, nch); eval(cmd);
    end

    %if n==size(files,1), save temp finalallspks; end
    
%     if isfield(spikes,'nsecs'),
%         finalallspks.nsecs = finalallspks.nsecs + spikes.nsecs;
%         cmd=sprintf('finalallspks.nsecs_%d = spikes.nsecs;',n);eval(cmd);
%     end
    
    if isfield(spikes,'nsecs_reduced'),
        finalallspks.nsecs_reduced = finalallspks.nsecs_reduced + spikes.nsecs_reduced;
        cmd=sprintf('finalallspks.nsecs_reduced_%d = spikes.nsecs_reduced;',n);eval(cmd);
        cmd=sprintf('finalallspks.nsecs_reduced_factor_%d = spikes.nsecs_reduced_factor;',n);eval(cmd);
    elseif isfield(spikes,'nsecs'),
        finalallspks.nsecs = finalallspks.nsecs + spikes.nsecs;
        cmd=sprintf('finalallspks.nsecs_%d = spikes.nsecs;',n);eval(cmd);
    end
  

    clear spikes

end % file loop

clear n files
clear ch
finalallspks.nsweeps=floor(max(finalallspks.fstimes)/1000);

spikes=finalallspks;
save ('joinedsortfile', 'spikes');


