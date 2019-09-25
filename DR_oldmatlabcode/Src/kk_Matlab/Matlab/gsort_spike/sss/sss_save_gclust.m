
function [savename] = sss_save_gclust(spikes, elec, day, date)
%%% CALL FROM GSORTM


gclustdata.params = spikes.matclustparams;
gclustdata.assigns = spikes.hierarchy.assigns;
clusters=unique(gclustdata.assigns);
clusters(find(clusters==0))=[];
gclustdata.clusters=clusters;

savename = ['gclust_' day '_' date '-' elec];
save(savename,'gclustdata');

