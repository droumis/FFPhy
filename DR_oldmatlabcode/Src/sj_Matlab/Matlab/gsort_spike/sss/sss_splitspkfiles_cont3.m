function [spikes,aspkfile]=sss_splitspkfiles_cont3(file,nspk,nsplitfiles)

%% sss_splitspkfiles_cont3('spiketetcut_all_thr-75',95000,2);
%% sss_splitspkfiles_cont3('spiketetcut_all_thr-34',109000,2);
%% Assuming 4 channels
ch=1:4;

load(file);
allspikes=spikes; clear spikes;
totalspks=length(allspikes.fstimes);

spkstart=1;
for n=1:nsplitfiles
    if n~=nsplitfiles,
        spkend=spkstart+(nspk-1);
    else
        spkend=totalspks;
    end
    
    spikes.waveforms=allspikes.waveforms(spkstart:spkend,:);
    for nch=ch(1):ch(4)
        cmd=sprintf('spikes.waveforms_ch%d = allspikes.waveforms_ch%d(spkstart:spkend,:);',nch,nch); eval(cmd);
    end
    spikes.spiketimes=allspikes.spiketimes(spkstart:spkend);
    spikes.swtimes=allspikes.swtimes(spkstart:spkend);
    spikes.fstimes=allspikes.fstimes(spkstart:spkend);
    spikes.ftimes=allspikes.ftimes(spkstart:spkend);
    
    spikes.nsweeps=allspikes.nsweeps;
    spikes.Fs=allspikes.Fs;
    spikes.threshT=allspikes.threshT;
    spikes.threshV=allspikes.threshV;
    
   
    
    filename=[file '-part' num2str(n)];
    save(filename,'spikes');
    
    spkstart=spkend+1;
end


