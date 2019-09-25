
function spikes = sss_cleanclustersfinal(spikes, nch)
%   CALL FROM GSORTM

%spktimes=spikes.fstimes;
neurons=spikes.hierarchy.assigns;
clustnames=unique(neurons);
ch=1:nch;

clustnames(find(clustnames==0))=[];
%%%%%% CLEAN %%%%%%%%%%%%%%
store_discard=[];

for clust = 1:length(clustnames)
    
    currclust=clustnames(clust)
    neurons=spikes.hierarchy.assigns;
   
    memberidxs = find(neurons==currclust);
    membertimes = spikes.fstimes(find(neurons==currclust));
    name=clustnames(clust);
    isis = diff(membertimes);
    tmin=1;      %ms
    tref=2;

    bin1 = find ((isis <= tmin) & (isis>0))+1; % number of coincident spikes
    bin2 = find ((isis <= tref) & (isis>tmin))+1;

%     d=randperm(length(bin2)); d=d(1:ceil(length(d)/3));
%     discard2 = memberidxs(bin2 (d)) ;
% 
%     d=randperm(length(bin1)); d=d(1:ceil(length(d)/1.5));
%     discard1 = memberidxs(bin1 (d)) ;
% 
%     discard = unique([discard1;discard2]);
%     store_discard=[discard];

    store_discard = unique([bin1; bin2]);

    spikes.fstimes(store_discard)=[];
    spikes.hierarchy.assigns(store_discard)=[];
    if isfield(spikes,'pca'), spikes.pca.scores(store_discard)=[]; end
    spikes.waveforms(store_discard,:)=[];
    spikes.spiketimes(store_discard)=[];
    spikes.swtimes(store_discard)=[];
    spikes.ftimes(store_discard)=[];
    for ncha=ch(1):ch(4)
        cmd=sprintf('spikes.waveforms_ch%d(store_discard,:) = [];',ncha); eval(cmd);
    end


    %         bin1 = (find ((isis <= tref) & (isis>0)));
    %         d=randperm(length(bin1)); d=d(1:floor(length(d)/2));
    %         discard1 = bin1(d) ;

end


