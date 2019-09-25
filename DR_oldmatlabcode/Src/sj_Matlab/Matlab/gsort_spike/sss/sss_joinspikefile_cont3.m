function [spikes,aspkfile]=sss_joinspikefile_cont3(spec1,spec2, ch, range_files, pushfiles)
%% sss_joinspikefile_cont3('spike_cut_','_tet_thr5', 4, 5:9, 0);
%% eg. spec1='spike_cut_'; spec2='_tet_thr5',ch=4; nfiles=10;
%% pushfiles=0 for no, 1 for yes

%%%% Inititalize allspikes

stimcode=0; unqstim=0;
if length(stimcode) < range_files(end)
    stimcode(1:range_files(end))=0;
end

allspikes.waveforms=[];
for nch=1:ch
    %cmd=sprintf('allspikes.waveforms_ch%d = [];',ch(nch)); eval(cmd);
    cmd=sprintf('allspikes.waveforms_ch%d = [];',nch); eval(cmd);
end
addsweep=0; addtrial=0; maxsweep=0;
allspikes.spiketimes=[]; allspikes.swtimes=[];
allspikes.fstimes=[];allspikes.ftimes=[];
allspikes.sweep = [];allspikes.trial = [];
allspikes.stimulus = []; allspikes.igorstim = [];allspikes.nsweeps = 0;
Alligorstim=[]; Allsweepidx=[]; sweepall=[];

%%%% Initialize allspikes

for n=range_files(1):range_files(end)
    n
    filename=[spec1 num2str(n) spec2]
    load(filename);

    if pushfiles==1
        spikes.spiketimes=spikes.ftimes/1000;
        spikes.fstimes=spikes.ftimes + 1000*(addsweep+1);
        spikes.swtimes=spikes.spiketimes+addsweep+1;
        spikes.sweep=floor((spikes.ftimes)/1000)+1+addsweep;
        maxsweep=max(spikes.sweep);
        spikes.trial=n*ones(size(spikes.ftimes));
        spikes.stimulus = stimcode(n) * ones(size(spikes.ftimes)); spikes.igorstim = stimcode(n) * ones(size(spikes.ftimes));
        newsweep = spikes.sweep + addsweep; newtrial = spikes.trial + addtrial;
        if autosweepstim==1, Allsweepidx=1:max(floor((spikes.ftimes)/1000)+1); sweepstim{n}=Allsweepidx+addsweep; end
        if autosweepall==1, sweepall=[sweepall Allsweepidx+addsweep]; end

        addsweep = maxsweep;
        addtrial = addtrial + 1;
    else
        Allsweepidx=1:max(floor((spikes.ftimes)/1000)+1); sweepstim{n}=Allsweepidx+addsweep;
        sweepall=[sweepall Allsweepidx+addsweep];
        maxsweep=max(spikes.sweep);  addsweep = maxsweep;
    end

    for nch=1:ch,
        cmd=sprintf('allspikes.waveforms_ch%d = [allspikes.waveforms_ch%d  spikes.waveforms_ch%d];',nch, nch, nch); eval(cmd);
    end

    allspikes.spiketimes=[allspikes.spiketimes;spikes.spiketimes]; allspikes.ftimes=[allspikes.ftimes;spikes.ftimes];
    allspikes.swtimes=[allspikes.swtimes;spikes.swtimes]; allspikes.fstimes=[allspikes.fstimes;spikes.fstimes];  %%% UPDATED FSTIMES
    allspikes.sweep = [allspikes.sweep; spikes.sweep]; allspikes.trial = [allspikes.trial; spikes.trial];
    allspikes.stimulus = [allspikes.stimulus; stimcode(n) * ones(size(spikes.ftimes))];
    allspikes.igorstim = [allspikes.igorstim; stimcode(n) * ones(size(spikes.ftimes))];
    allspikes.nsweeps = length(range_files);

    Fs=spikes.Fs; spikestarts=spikes.threshT;
    thresh=spikes.threshV(1);
    if thresh~=-Inf, thrs=thresh; else, thrs=spikes.threshV(2); end
    stimonset=spikes.stimonset; bckwindow=spikes.window; sweepd=spikes.sweepd;

    clear spikes
    
end % file loop

save temp
if (exist('sweepall')), allspikes.sweepall=sweepall; end
if (exist('sweepstim'))
    for i=1:length(unqstim), cmd=sprintf('allspikes.sweep_%d = sweepstim{i};',unqstim(i)); eval(cmd); end
end
if (exist('Alligorstim')), allspikes.Alligorstim=Alligorstim; end

for nch=1:ch
    cmd=sprintf('allspikes.waveforms = [ allspikes.waveforms; allspikes.waveforms_ch%d];',nch); eval(cmd);
end
allspikes.waveforms = allspikes.waveforms;

spikes=allspikes;

spikes.Fs=Fs; spikes.threshT=spikestarts;
if thrs<0, spikes.threshV=[thrs Inf]; else, spikes.threshV=[-Inf thrs];end
spikes.stimonset=stimonset; spikes.window=bckwindow;spikes.sweepd=sweepd;
spikes.nsweeps=floor(max(spikes.fstimes)/1000);

aspkfile = ([spec1 spec2 '_all']);
save (aspkfile, 'spikes');
