
function [indices] = sss_savespksforMclust(spikes,nspks, nch)
%%% CALL FROM GSORTM
%% TO only keep subset of spikes for statistics in Mclust
% sss_savespksforMclust('usesortf1_rename',20000);

ch=1:nch;

spikes.waveforms(nspks+1:end,:)=[];
spikes.spiketimes(nspks+1:end)=[];
spikes.swtimes(nspks+1:end)=[];
spikes.fstimes(nspks+1:end)=[];
spikes.ftimes(nspks+1:end)=[];
for nch=ch(1):ch(4)
    cmd=sprintf('spikes.waveforms_ch%d(nspks+1:end,:) = [];',nch); eval(cmd);
end

spikes.hierarchy.assigns(nspks+1:end)=[];
if isfield(spikes,'pca'),
    spikes.pca.scores(nspks+1:end)=[];
end

assigns=spikes.hierarchy.assigns;
nclusters = length(unique(assigns));

assigns=[nclusters;assigns];  %%% Mclust loads Klustakiwk from 2nd elemnt onwards
save('TT_tetNF_gsort.clu.1','assigns','-ascii');

%spikes=rmfield(spikes,{'overcluster';'sweep';'trial';'stimulus';'igorstim';'nsweeps';'outliers'});

if isfield(spikes,'pca'),
    spikes=rmfield(spikes,{'pca'});
end

savename= ['TT_tetNF_gsort_for_mclust'];
save (savename, 'spikes');

