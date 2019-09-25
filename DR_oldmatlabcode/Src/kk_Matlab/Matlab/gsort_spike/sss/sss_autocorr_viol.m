
function spikes = sss_isi_viol(spikes, useassigns, show, plot_ampl);
% figure; hold on; sss_finalplot(spikes, spikes.hierarchy.assigns, [1,2,15], 0.25);
%    OR CALL FROM Gsortm

if (nargin < 2)
    useassigns = spikes.hierarchy.assigns;
end
if (nargin < 3)
    show = unique(useassigns);
end
if (nargin < 4)
    plot_ampl = -0.3;
end

show = reshape(show, 1, []);
show(find(show==0))=[];

clf;
ylims = [min(spikes.waveforms(:)) max(spikes.waveforms(:))];


%set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold'); set(0,'defaultaxeslinewidth',2);
for clust = 1:length(show)
    
    members = find(useassigns == show(clust));
    memberwaves = spikes.waveforms(members,:);
    membertimes = sort(spikes.fstimes(members));
    
    if (size (memberwaves,1) > 1000)
        unit_sm = memberwaves(randperm (size (memberwaves,1)),:);
        unit_sm = unit_sm (1:1000, :);
        %plot(unit_sm',clr(1)); hold on;
        %plot (mean (tempwaves), clr(4), 'Linewidth', 2);
    else
        unit_sm=memberwaves;
    end

    tmin=1;
    tref = max(2, tmin*1.5);
    %isis = sort(diff(membertimes));
    isis = diff(membertimes);
    nrcsp = length(find ((isis <= tmin) & (isis>0))); % number of coincident spikes
    viol = find ((isis <= tmin) & (isis>0)); % idxs coincident spikes
    
    cnt=0; remove_idxs=[];
    for i=1:length(viol),
        
        f1=figure; redimscreen70s; 
        hold on;
        plot(unit_sm','b'); hold on;
        plot (mean (memberwaves),'y','Linewidth', 4);
        plot(memberwaves(viol(i),:),'r','Linewidth',2);
        if viol>1
            plot(memberwaves(viol(i)+1,:),'g','Linewidth',2);
        end
        
        %Choose waht to remove: red=viol=1, green=viol-1=2, none=0, both =3
        Tlines={['Spike: 0=None, 1=Red, 2=green, 3=both']};
        spkidx = inputdlg(Tlines, 'Enter Spike idx', 1, {'1'});
        if ~isempty(spkidx)              %see if user entered something and did not hit cancel
            stimn = spkidx{1};
            stimn=str2num(stimn);
        else
            stimn=0;
        end
        
        switch stimn
            case 1
                cnt=cnt+1;
                remove_idxs(cnt)=viol(i);
            case 2
                cnt=cnt+1;
                remove_idxs(cnt)=viol(i)+1;
            case 3
                cnt=cnt+1;
                remove_idxs(cnt)=viol(i)+1;
                remove_idxs(cnt)=viol(i);
            case 0
                cnt=cnt;
        end
        
        close(f1);
    end
    
   
    
   
    
    %%% Update assignments structure by pushing removed spikes to 0
    useassigns(members(remove_idxs)) = 0;
    spikes.hierarchy.assigns = useassigns;
    %membertimes(remove_idxs)=[];
    %memberwaves(remove_idxs)=[];
    
    
    
    %Get the cluster again and plot new auto-correlation
    
    members = find(useassigns == show(clust));
    memberwaves = spikes.waveforms(members,:);
    membertimes = sort(spikes.fstimes(members));
    isis = diff(membertimes);
    nrcsp = length(find ((isis <= tmin) & (isis>0))); % number of coincident spikes
    
    figure; hold on;
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
    
    xlabel(['ms (' num2str(acorr_bin_msec) 'ms binsize)']);
    ylabel('Rate');
    
    %if clust==1
        msgstr = {};
        msgstr{end+1} = 'New AUTOCORR ';
       	msgstr{end+1} = ['Coinspks= ' num2str(nrcsp) '; %= ' num2str( roundn( 100*nrcsp/length(membertimes),-2))];   
        title(msgstr)
    %else
    %    title(['Coinspks= ' num2str(nrcsp) '; %= ' num2str( roundn( 100*nrcsp/length(membertimes),-2))] );
    %end
     drawnow;

    
    

end
