
%   'window'   -- 1x2 vector specifies the window before and after each included event.
%                   Default is 500 mseconds before and after
%                   event start time.
%   'binsize'  -- histc binsize for spike rasters, 1 ms


function out = dfa_licktrigspiking(index, excludeperiods, varargin)

fprintf('%d %d %d %d\n',index)

eventName = '';
window = [0.5 0.5]; % in sec
binsize = 0.001; % 1 ms for rasters
plotfigs = 0;
TF = 1; % tetfilter number (legacy)
events = [];
if ~isempty(varargin)
    assign(varargin{:});
end

out.index = index;
out.psth = [];
out.dims = {'eventnum', 'time'};
out.eventtimes = [];
out.time = [];
out.bin = binsize;
out.window = window;

day = index(1);
epoch = index(2);
ntrode = index(3);
cluster = index(4);

try
    spiketimes = spikes{day}{epoch}{ntrode}{cluster}.data(:,1);
catch
    return
end

% get events
if isempty(events)
    if any(eventName)
        evvar = varargin{find(cellfun(@(x) ~isempty(x), ...
            strfind(varargin(1:2:end), eventName), 'un', 1))*2-1};
    else
        evvar = varargin{find(cellfun(@(x) ~isempty(x), ...
            strfind(varargin(1:2:end), 'kons'), 'un', 1))*2-1};
    end
    events = eval(evvar);
end
% First receive valid event periods for the epoch's day.
try
    ec = events{day}{epoch};
catch
    ec = events{day}{epoch}{TF};
end
eventtimes = ec.starttime;
% eptimes = ec.timerange(1):1/ec.samprate:ec.timerange(end);

% includetimes
includetimes = ~isExcluded(eventtimes, excludeperiods); % list of ones and zeros sampled every millisecond, ones = included, zeros = excluded
eventtimes = eventtimes(includetimes);
if isempty(eventtimes)
    return
end
% includetimes = includetimes(:);
% if sum(includetimes)==0
%     return
% end

% histc bins, where center bin is centered time 0
time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
if ~isempty(eventtimes)
    psth = nan(length(eventtimes),length(time));
    for r=1:length(eventtimes)
        shist = histc(spiketimes, eventtimes(r) + time);
        psth(r,:) = shist;
    end
end

out.index = index;
out.psth = psth;
out.dims = {'eventnum', 'time'};
out.eventtimes = eventtimes;
out.time = time;
out.bin = binsize;
out.window = window;

if plotfigs
    [xx, yy] = find(out.psth');
    % f1 = scatter(xx/1000-1.001,yy,10, '+k');
    % f1.MarkerEdgeAlpha = .4;
    figure
    s1 = scatterhist(xx/1000-1.001,yy, 'Kernel', 'on', 'bandwidth', ...
        [.015; length(eventtimes)*.01], 'location','SouthEast', 'Direction', 'in', 'Marker', ...
        '+', 'color', 'k', 'MarkerSize', 2);
    s1(3).Position(3) = .08;
    s1(2).Position(4) = .08;
    ylabel('licknum','FontSize',8,'FontWeight','bold', 'FontName','Arial')
    xlim([-window(1) window(2)]); xticks([-window:.2:window]);
    xlabel('time s','FontSize',8,'FontWeight','bold', 'FontName','Arial');
    axis tight
    line([0 0],ylim, 'color','red', 'linewidth', 1)
end
end




