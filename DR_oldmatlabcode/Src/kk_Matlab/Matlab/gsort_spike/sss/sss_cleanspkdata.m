
function spikes = sss_cleanspkdata(filename,removeclu);
%% TO remove specific Clusters from spike data and refresh it for another round of the algorithm

%sss_cleanspkdata('spiketetcut_all_thr-25_sort2.mat',[0,32]);
%sss_cleanspkdata('spiketetcut_all_thr-30_sort1',[0,19,32]);

load(filename);
ch=1:4;

% noiseidx=spikes.noiseidx;
% spikes.waveforms(noiseidx,:)=[];
% spikes.spiketimes(noiseidx,:)=[];
% spikes.swtimes(noiseidx,:)=[];
% spikes.fstimes(noiseidx,:)=[];
% spikes.ftimes(noiseidx,:)=[];
% for nch=ch(1):ch(4)
%     cmd=sprintf('spikes.waveforms_ch%d(noiseidx,:) = [];',nch); eval(cmd);
% end

noiseidx=[];
for i=1:length(removeclu)
    noiseidx = [noiseidx; find(spikes.hierarchy.assigns==removeclu(i))];
end

spikes.waveforms(noiseidx,:)=[];
spikes.spiketimes(noiseidx,:)=[];
spikes.swtimes(noiseidx,:)=[];
spikes.fstimes(noiseidx,:)=[];
spikes.ftimes(noiseidx,:)=[];
for nch=ch(1):ch(4)
    cmd=sprintf('spikes.waveforms_ch%d(noiseidx,:) = [];',nch); eval(cmd);
end

spikes = rmfield(spikes, {'noiseidx';'tictoc';'overcluster';'hierarchy';'parameters'});

savename= [filename '_clean'];
save (savename, 'spikes');

