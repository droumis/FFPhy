function [spikes,aspkfile]=sss_joinspikefile_cont5(spec1,spec2, ch, range_files, pushfiles)
%% sss_joinspikefile_cont5('spike_cut_','_tet_thr-34', 4, 1:8, 0);
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
%allspikes.sweep = [];allspikes.trial = [];
%allspikes.stimulus = []; allspikes.igorstim = [];
allspikes.nsweeps = 0;
%Alligorstim=[]; Allsweepidx=[]; sweepall=[];

%%%% Initialize allspikes

for n=range_files(1):range_files(end)
    n
    filename=[spec1 num2str(n) spec2]
    load(filename);

    if pushfiles==1
        spikes.spiketimes=spikes.ftimes/1000;
        spikes.fstimes=spikes.ftimes + 1000*(addsweep);
        spikes.swtimes = spikes.fstimes/1000;
        spikes.sweep=floor((spikes.ftimes)/1000)+addsweep;
        maxsweep=max(spikes.sweep);
   
        %addsweep = maxsweep;
         if size(data,1)<Fs*300      %%% LAST FILE ADDSWEEP WILL GIVE THE SECOND FOR LAST SPIKE TIME 
                addsweep = maxsweep+1;
         else
                addsweep = addsweep + 300;
         end               
        addtrial = addtrial + 1;
    else
    end

%     for nch=1:ch,
%         cmd=sprintf('allspikes.waveforms_ch%d = [allspikes.waveforms_ch%d  spikes.waveforms_ch%d];',nch, nch, nch); eval(cmd);
%     end

    fullwave=[];
    for nch=1:ch,
        cmd=sprintf('fullwave = [fullwave;  spikes.waveforms_ch%d];',nch); eval(cmd);
    end
    
    allspikes.waveforms = [ allspikes.waveforms fullwave];   
    
    allspikes.spiketimes=[allspikes.spiketimes;spikes.spiketimes]; allspikes.ftimes=[allspikes.ftimes;spikes.ftimes];
    allspikes.swtimes=[allspikes.swtimes;spikes.swtimes]; allspikes.fstimes=[allspikes.fstimes;spikes.fstimes];  %%% UPDATED FSTIMES
  
    Fs=spikes.Fs; spikestarts=spikes.threshT;
    thresh=spikes.threshV(1);
    if thresh~=-Inf, thrs=thresh; else, thrs=spikes.threshV(2); end
    %stimonset=spikes.stimonset; bckwindow=spikes.window; sweepd=spikes.sweepd;

    clear spikes
    
end % file loop





spikes=allspikes;

spikes.Fs=Fs; spikes.threshT=spikestarts;
if thrs<0, spikes.threshV=[thrs Inf]; else, spikes.threshV=[-Inf thrs];end

aspkfile = (['spiketetcut_all' spec2 ]);
save (aspkfile, 'spikes');
