function plotraster5(timevec,psth,markersize,color,varargin)

%% Note that this plot FLATTENS any bins with more than one spike. Treats any nonzero bins a single spike.
%% This is actually quite important if dfskk_runriptriggeredspiking2 has a binsize bigger than 2 ms (can contain two spikes)

%plots raster given a certain format of data :

% timevec is a histc-style vector of times
% psth is a matrix of 0s and 1s, with each row vector corresponding to a
        % trial
% height is the height of the raster lines
% if burstisi is specified, will highlight spikes that have an isi, before or after,
    %   within burstisi (in ms) 

burstisi = 0;
postripmatrix = [];

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'burstisi'               % the burst threshold for which to highlight bursts (bigger dot)
            burstisi = varargin{option+1};
        case 'postripplematrix'
            postripplematrix = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

timevec = timevec*1000;  % convert to ms
% obtain bincenter form
bincenters = timevec(1:(end-1)) + (timevec(2)-timevec(1)) / 2 ;
% also truncate last bin of histc psth output
psth = psth(:,1:(end-1));


% burst psth %%%%%%%%%%%%%%%%%%%%%
    % if specified, plot burst-associated spikes
    if 0
    if burstisi > 0
        
        % initialize output
        burstspike_psth = zeros(size(psth));  % notice we're using the original psth to detect bursts
        
        % iterate through each trial
        for tt = 1:size(psth,1)
            if sum(psth(tt,:)) > 1
                % first identify isis of spikes in this trial
                spikepositions = find(psth(tt,:));
                prespike_isis = diff([-inf  spikepositions  inf]);
                prespike_isis = prespike_isis(1:(end-1));
                postspike_isis = [prespike_isis(2:end) inf];
                validinds = (prespike_isis <= burstisi) | (postspike_isis <= burstisi);
                burstspike_psth(tt,spikepositions(validinds)) = 1;
            end
        end
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot    %%%%%%%%%%%%%%%%%%%

for spikeset = [ 1 2 1]
    
    if spikeset == 1       % all spikes
        clr = color;
        plotpsth = psth;
    elseif spikeset == 2       % ripple-associated spikes
        if ~isempty(postripplematrix) 
            ripmatrix = [ zeros(size(postripplematrix))  ones(size(postripplematrix,1),1)  postripplematrix ];
            plotpsth = psth.*ripmatrix;
            ypoints2 = ripmatrix.*levelmat;   
            clr = [.5 .5 .5];
            %clr = [.1 .1 .1]
            %clr = color;
        end
    end
    
    xpoints = plotpsth.*repmat(bincenters,size(plotpsth,1),1);
        xpoints(xpoints==0)=[];
    xpoints = xpoints(:);
    levelvec = (1:size(plotpsth,1))';
    levelmat = repmat(levelvec,1,length(bincenters));
    ypoints = plotpsth.*levelmat;
        ypoints(ypoints==0)=[];
        ypoints = ypoints(:);
    
    if spikeset == 1  % regular spike plot
        scatter(xpoints,ypoints,markersize,clr,'.')
    elseif spikeset == 2  % ripple-associated spike
        if 0
        hold on
        scatter(xpoints,ypoints,markersize+100,[1 .65 .65],'.')    % highlight ripple-associated spikes
        end
    end
    
    % burst
    if 0
    if burstisi > 0
        % filter for either all spikes or rip-spikes
        plotburstpsth = burstspike_psth & plotpsth;
            xpoints = plotburstpsth.*repmat(bincenters,size(plotburstpsth,1),1);
                xpoints(xpoints==0)=[];
                xpoints = xpoints(:);
            levelvec = (1:size(plotburstpsth,1))';
            levelmat = repmat(levelvec,1,length(bincenters));
            ypoints = plotburstpsth.*levelmat;
                ypoints(ypoints==0)=[];
                ypoints = ypoints(:);
        hold on
        %burst spike plot
        scatter(xpoints,ypoints,markersize+150,clr,'.')
    end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
% axes
xlim([bincenters(1) bincenters(end)])
ylim([0 size(psth,1)])
set(gca,'YDir','reverse')
    % x ticks
%set(gca,'XTick',mean([timevec(1) timevec(2)]):500:mean([timevec(end-1) timevec(end)]))
set(gca,'XTick',[])
% label axes
set(gca,'FontSize',18,'FontWeight','bold')
%xlabel('time (ms)','FontSize',16,'FontWeight','bold')
ylabel('ripple #','FontSize',18,'FontWeight','bold')

end