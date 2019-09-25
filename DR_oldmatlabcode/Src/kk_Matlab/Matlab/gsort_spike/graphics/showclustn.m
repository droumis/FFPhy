
function sss_finalplot(spikes, useassigns, show, plot_ampl);
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
    
    subplot(length(show),3, 3 * (clust-1) + 1);
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
    
    subplot(length(show),3,3 * (clust-1) + 2);
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
    
    %% SJ: Apr 16
    membertimes_ms = sort(spikes.fstimes(members)); isis_ms = sort(diff(membertimes_ms)); isis_ms = isis_ms(isis_ms < 0.01);
    
   if length(isis)==0, isis=0.009; end
   n = histc (isis,0:0.001:0.01); 
   bar(0:0.001:0.01,n,'histc');
%     [n,x] = hist(isis,10);
%     plot(x,n,'.-');
    set(gca,'Xlim',[0 0.01]);
    if (clust < length(show))
        set(gca,'XTickLabel',{});
    else
        xlabel('ISI (sec)');
    end
    title(['ISI score = ' num2str(scores(1)) '; CoinSpks= ' num2str(nrcsp)]);
    clear isis n
    
    subplot(length(show),3,3 * (clust-1) + 3);
    isis = sort(diff(membertimes));
    isis = isis(isis < 2);
    if length(isis)==0, isis=1.9; end
    n = histc (isis,0:0.02:2); 
    bar(0:0.02:2,n,'histc');
%     [n,x] = hist(isis,100);
%     plot(x,n,'.-');
    set(gca,'Xlim',[0,2]);
    if (clust < length(show))
        set(gca,'XTickLabel',{});
    else
        xlabel('ISI(sec)');
    end
    title('ISI Histogram');
end