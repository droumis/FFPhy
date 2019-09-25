
function sss_densityplot(spikes, useassigns, show, plot_ampl)
%    CALL FROM Gsortm

if (nargin < 2)
	useassigns = spikes.hierarchy.assigns;
end
if (nargin < 3)
    show = unique(useassigns);
end
if (nargin < 4)
    plot_ampl = 300;
end

show = reshape(show, 1, []);

clf;
ylims = [min(spikes.waveforms(:)) max(spikes.waveforms(:))];
%figure; hold on; 
redimscreen 
%redimscreen_halfvert(0);
a=50; b=200; c=1000;

for clust = 1:length(show)
    members = find(useassigns == show(clust));
    memberwaves = spikes.waveforms(members,:);
 
    membertimes = sort(spikes.swtimes(members));
    % MY CHANGE
%     membertimes = (spikes.spiketimes(members));
    
    subplot(length(show),4, 4 * (clust-1) + 1);
    [n,x,y] = hist2d(memberwaves);
    imagesc(x,y,n); axis xy; colormap hot;
    if (clust < length(show))
        set(gca,'XTickLabel',{});
    end
    %set(gca, 'YLim', ylims);
    set(gca, 'YLim', [-plot_ampl plot_ampl]);
    title(['# Spikes = ' num2str(size(members,1))]);
    if (show(clust) ~= 0)
        ylabel(['Cluster# ' num2str(show(clust))]);
    else
        ylabel('Outliers');
    end
    
    subplot(length(show),4,4 * (clust-1) + 2);
	tmin = size(spikes.waveforms,2)./spikes.Fs;
	tref = max(0.002, tmin*1.5);
    tmin=0.001;      % MY CHANGE
    tref=0.002;     % MY CHANGE
    [a, scores] = isiQuality(membertimes, membertimes, tmin, 0.010, tref, spikes.Fs);
    isis = sort(diff(membertimes));
    
    %nrcsp = length(find(isis <= tmin)); % number of coincident spikes 
%     x=find ((isis <= tmin));
%     isis(x(1:end-2))=[];
    nrcsp = length(find ((isis <= tmin) & (isis>0))); % number of coincident spikes 
    
    isis = isis(isis < 0.01);
%   isis = isis(find ((isis < 0.01) & (isis>0)));  % MY CHANGE
    
%     %% SJ: Apr 16
    %%% Jan09
    membertimes_ms = sort(spikes.fstimes(members)); 
    isis_ms = sort(diff(membertimes_ms));
     
    
   %ISIS-A 
   isis_a = isis_ms(find(isis_ms < 20));  %%% 20 ms  
   if length(isis_a)==0, isis_a=19; end
   n = histc(isis_a,0:1:20);
   %n(1)=0;
   norm_n = n./max(n);
   bar([0:1:20],n,'histc');
%     [n,x] = hist(isis,10);
%     plot(x,n,'.-');
    set(gca,'Xlim',[-2 20]);
    %set(gca,'Ylim',[0 1.05]);
   % if (clust < length(show))
   %     set(gca,'XTickLabel',{});
   % else
        xlabel('ISI (ms)');
   % end
    title(['ISI score = ' num2str(scores(1)) '; CoinSpks= ' num2str(nrcsp)]);
    %clear isis n
    hold on
    plot([1 1], get(gca, 'YLim'), 'r:','LineWidth',2)
    hold off
    
    
    
    
    %ISIS-B
    subplot(length(show),4,4 * (clust-1) + 3);
    isis_b = isis_ms(find(isis_ms < 10000));  %%% Long Time Scale in  ms  
    if length(isis_b)==0, isis_b=9999; end

    n = histc(isis_b,0:10:10000); 
    %n([1,2])=0;
    norm_n = n./max(n);
    plot([0:10:10000], n,'k-','Linewidth',2);    
    set(gca,'Xlim',[-10,1000]);
    %set(gca,'Ylim',[0 1.05]);
    %set(gca, 'XScale', 'log', 'XLim', [10^(-1) 10^4]);
    %set(gca, 'XScale', 'log', 'XTick', [10^(0),10^(2),10^(4)]);
    hold on
    plot([1 1], get(gca, 'YLim'), 'r:','LineWidth',2)
    hold off
    
    %if (clust < length(show))
    %    set(gca,'XTickLabel',{});
    %else
        xlabel('ISI(ms)');
    %end
    title('ISI Histogram');
    
    
    
    
    %%%% HISTOGRAM OF LOG(ISI) FROM MCLUST
   
    tMclust = membertimes_ms*10;  % timestamps in 0.1ms resolution
    tsMcl = ts(tMclust,ts);
    [H, binsUsed] = HistISI(tsMcl);
    normH = H./max(H);
    subplot(length(show),4,4 * (clust-1) + 4);
    plot(binsUsed, normH);
    set(gca, 'XScale', 'log', 'XLim', [10^(-1) 10^3]);
    set(gca, 'XScale', 'log', 'XTick', [10^(0),10^(2),10^(3)]);
    set(gca, 'Ylim', [0 1.2*max(normH)]);
     title('histogram of log(ISI)');
    ylabel('nSpikes');
    xlabel('msec');
     hold on
    plot([1 1], get(gca, 'YLim'), 'r:','LineWidth',2)
    hold off
    

end

%keyboard;
