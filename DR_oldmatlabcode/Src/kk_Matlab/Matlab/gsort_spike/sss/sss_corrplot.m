
function sss_corrplot(spikes, useassigns, show, plot_ampl, nch)
% figure; hold on; sss_finalplot(spikes, spikes.hierarchy.assigns, [1,2,15], 0.25);
%    OR CALL FROM Gsortm

if (nargin < 2)
    useassigns = spikes.hierarchy.assigns;
end
if (nargin < 3)
    show = unique(useassigns);
end
if (nargin < 4)
    plot_ampl = 300;
end
if (nargin < 5)
    nch = 4;
end

if isfield(spikes,'nsecs_reduced'),
    nsecs=spikes.nsecs_reduced;
elseif isfield(spikes,'nsecs'),
    nsecs=spikes.nsecs;
end

show = reshape(show, 1, []);
show(find(show==0))=[];

clf;
ylims = [min(spikes.waveforms(:)) max(spikes.waveforms(:))];
figure; hold on;
redimscreen;
%redimscreen_halfvert(0);

set(0,'defaultaxesfontsize',12);
set(0,'defaultaxesfontweight','bold'); set(0,'defaultaxeslinewidth',2);

totalclu=length(show);
for clust = 1:length(show)
    
    % Current Cluster 
    members = find(useassigns == show(clust));
    memberwaves = spikes.waveforms(members,:);
    membertimes = sort(spikes.fstimes(members));
    
    tmin=1;
    tref = max(2, tmin*1.5);
    %isis = sort(diff(membertimes));
    isis = diff(membertimes);
    nrcsp = length(find ((isis <= tmin) & (isis>0))); % number of coincident spikes
    
    % Autocorr for current cluster
    subplot(totalclu,totalclu,(clust-1)*length(show)+clust); hold on;
    
    tMclust = membertimes*10;  % timestamps in 0.1ms resolution
    tsMcl = ts(tMclust,ts);
    acorr_bin_msec = 1; scale=30; %ms
    [histvals,x] = AutoCorr(Data(tsMcl), acorr_bin_msec, scale); % This is peter's C AutoCorr
    
    px=[flipud(-x);x]; phistvals=[flipud(histvals); histvals];
    bar(px,phistvals,'FaceColor','k');			% show acorr
    m=1.2*max(histvals);  arr=0:0.1:m;
    plot(1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
    plot(-1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
        axis([-scale scale 0 1.3*m]);
    if (clust == 1)
        xlabel(['ms (' num2str(acorr_bin_msec) 'ms binsize)'],'FontSize',14,'Fontweight','bold');
        ylabel('Rate','FontSize',14,'Fontweight','bold');
    end
    
    msgstr = {};
    msgstr{end+1} = ['Clu ' num2str(show(clust))];
    msgstr{end+1} = ['Coinspks= ' num2str(nrcsp) '; %= ' num2str( roundn( 100*nrcsp/length(membertimes),-2))] ;   
    title(msgstr,'FontSize',14,'Fontweight','bold');
    
    drawnow;
      
    
    for i=clust+1:length(show),    
        
        cluidx = find(useassigns == show(i));
        clutimes = sort(spikes.fstimes(cluidx));
        
        % Clust vs i Cross-Corr
        subplot(totalclu,totalclu,(clust-1)*length(show)+i); hold on;
        
        clu_tsMcl = ts(clutimes*10,ts);
        [chistvals,cx] = CrossCorr(Data(tsMcl),Data(clu_tsMcl),acorr_bin_msec,2*scale); %[C, B] = CrossCorr(t1,t2,binsize,nbins);
        
        bar(cx,chistvals,'FaceColor','k');			% show acorr
        m=1.2*max(chistvals);  arr=0:0.1:m;
        plot(1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
        plot(-1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
        axis([-scale scale 0 1.3*m]);
        
        title([num2str(show(clust)) ' vs ' num2str(show(i))],'FontSize',12,'Fontweight','bold')
        
        %crossisi = membertimes-clutimes;
        %nrcsp = length(find ((isis <= 1) & (isis >= -1)));
        %title(['CrCoinspks= ' num2str(nrcsp) '; %= ' num2str( roundn( 100*nrcsp/length(membertimes),-2))] );
        drawnow;
    
        
        % i vs. Clust Cross-Corr
        
        subplot(totalclu,totalclu,(i-1)*length(show)+clust); hold on;
        
        clu_tsMcl = ts(clutimes*10,ts);
        [chistvals,cx] = CrossCorr(Data(clu_tsMcl),Data(tsMcl),acorr_bin_msec,2*scale); %[C, B] = CrossCorr(t1,t2,binsize,nbins);
        
        bar(cx,chistvals,'FaceColor','k');			% show acorr
        m=1.2*max(chistvals);  arr=0:0.1:m;
        plot(1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
        plot(-1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
        axis([-scale scale 0 1.3*m]);
        
        title([num2str(show(i)) ' vs ' num2str(show(clust))],'FontSize',12,'Fontweight','bold');
        drawnow;
        
               
    end
    
end
