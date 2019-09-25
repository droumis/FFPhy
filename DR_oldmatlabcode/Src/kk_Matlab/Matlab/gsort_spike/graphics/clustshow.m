function clustshow(spikes, useassigns, show,nch)
%    CLUSTSHOW  another temporary utility to show clusters
%       CLUSTSHOW(SPIKES, [USEASSIGNS], [SHOW]);

if (nargin < 2)
	useassigns = spikes.overcluster.assigns;
end

if ((nargin < 3) || isempty(show))
    show = unique(useassigns);
end

if (nargin < 4)
	nch=4;
end


if (length(show) < 10)
	cmap = winter(256);
else
    cmap = jet(256);
end
cind = floor(linspace(1,256,length(show)));
cind = cind(randperm(length(cind)));

%%%%%%%%%%%%% Height & FWHM estimation
if nch ==1
    [heights,widths] = thresholded_peaks(spikes);
else
    [heights,widths] = thresholded_peaks_fortet(spikes,spikes.waveforms);
end
widths = widths + 0.3 * randn(size(widths));   % jitter the widths a little bit

%%%%%%%%%%%%%%%

subplot(2,1,1); cla; hold on; set(gca,'xlim',[1 32]);
subplot(2,1,2); cla; hold on; legend off;
for k = 1:length(show)
    members = find(ismember(useassigns,show(k)));
    if (show(k) ~= 0)
        symbol = '.';
    else
        symbol = 'x';
    end
    if (~isempty(members))
        subplot(2,1,1);
        h = mplot(spikes.waveforms(members,:), 'Color', cmap(cind(k),:));
        set(h, 'ButtonDownFcn', {@raise_me, h});
        subplot(2,1,2);
        h = plot(heights(members),widths(members),symbol,'Color',cmap(cind(k),:));
        set(h, 'ButtonDownFcn', {@raise_me, h});
    end
end

subplot(2,1,1); hold off; axis tight; xlabel('Time (samples)');  ylabel('Voltage (A/D Levels)');
subplot(2,1,2); hold off; xlabel('Height (A/D Levels)');  ylabel('Width (samples)');
if (length(show) < 8)
    leg = cell(length(show),1);
    for k = 1:length(show)
        leg{k} = num2str(sort(show(k))); 
    end
    legend(leg,0);        
end
