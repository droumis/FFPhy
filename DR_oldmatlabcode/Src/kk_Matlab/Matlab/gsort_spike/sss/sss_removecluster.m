
function [spikes] = chr_removecluster(spikes,cluname)

%%% In Gsortm : To Remove spikes of a cluster in dataset
%% TO only keep subset of spikes for statistics in Mclust
% chr_reducespksforMclust('spiketetcut_all_thr-100-red-sortf1_rename.mat',20000)

ch=1:4;

idxs=find(spikes.hierarchy.assigns==cluname);

spikes.waveforms(idxs,:)=[];
spikes.spiketimes(idxs)=[];
spikes.swtimes(idxs)=[];
spikes.fstimes(idxs)=[];
spikes.ftimes(idxs)=[];
for ncha=ch(1):ch(4)
    cmd=sprintf('spikes.waveforms_ch%d(idxs,:) = [];',ncha); eval(cmd);
end

spikes.hierarchy.assigns(idxs)=[];
if isfield(spikes,'pca'), spikes.pca.scores(idxs)=[]; end


