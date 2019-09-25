function sss_drift(spikes, assignments, show, nch, usech)
% Shantanu - Apr 2011: to plot drift in amplitude for clusters

if (isempty(assignments))
    if (isfield(spikes, 'hierarchy') && isfield(spikes.hierarchy, 'assigns'))
        assignments = spikes.hierarchy.assigns;
    elseif (isfield(spikes, 'overcluster'))
        assignments = spikes.overcluster.assigns;
    else
        assignments = ones(size(spikes.spiketimes));
    end
end

if (isempty(show)),  show = clusters;  end;
if (nargin < 4), nch=1; end
if (nargin < 5), usech=1; end

%clusters = unique(assignments); 
clusters = show; numclusts = length(clusters); 
%cmap=jet(length(show));
cmap = hot(2*length(show));


channel_lth = size(spikes.waveforms,2)/nch;
nSpikes = size(spikes.waveforms,1);

%%% GETTING AMPL FOR DRIFT - Use only first channel if multiple channels %%
if nch==1,
    ampl=abs(max(spikes.waveforms,[],2)) + abs(min(spikes.waveforms,[],2));
else
    % Calculate Amplitude for mulit-channel data
    for ch=1:nch
        if isfield(spikes,'waveforms_ch1')
            cmd=sprintf('ampl(:,ch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',ch,ch); eval(cmd);
        else
            w = spikes.waveforms(:,channel_lth*(ch-1)+1:channel_lth*ch);
            ampl(:,ch) = abs(max(w,[],2)) + abs(min(w,[],2));
        end
    end
end


%%% TO NOT SHOW OUTLIERS
show(show == 0) = [];
data=ampl(:,usech);
spktimes = spikes.fstimes./10000;

% Make the plot.
figure;  hfig = gcf;  hax = gca;    hold on;
%redimscreen70s;
%orient(gcf,'landscape');
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    
for clu = 1:numclusts
    sel = find(assignments == clusters(clu));
    inds{clu} = sel;
    hndl = plot(spktimes(sel), data(sel), '.','MarkerSize',6);
    set(hndl, 'ButtonDownFcn', {@raise_me, hndl});

    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.1);
    end;
	
	if (clusters(clu) == 0),                 h_out = hndl;
	elseif (ismember(clusters(clu), show)),  h_pts(clusters(clu)) = hndl;
	end
	
	ssghandle(clu) = hndl;
end

h_pts = h_pts(show);
if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;

ylabel( ['Ampl (uV) on Ch' num2str(usech)]); xlabel('Time (sec)');
hold off;

