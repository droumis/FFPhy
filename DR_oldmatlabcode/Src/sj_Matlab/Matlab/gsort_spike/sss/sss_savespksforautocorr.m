
function sss_savespksforautocorr(spikes, saveclusters)
%%% CAN CALL FROM GSORTM

if nargin<2,
    saveclusters=unique(spikes.hierarchy.assigns);
end

assigns=spikes.hierarchy.assigns;
clusters=saveclusters;
clusters(find(clusters==0))=[];

for i=1:length(clusters)
    spkdata.spktimes{i}=spikes.fstimes(find(assigns==clusters(i)));
    spkdata.neurons{i}=clusters(i);
end

if isfield(spikes,'nsecs')  
    spkdata.nsecs=spikes.nsecs;
end

savename= ['Sort_autocorr'];
save (savename, 'spkdata');

