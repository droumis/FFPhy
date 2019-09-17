



function out = dfa_lickswrcorr(idx, excludeIntervals, varargin)
pconf = paramconfig;
displayplots = 0;
saveplots = 1;
savefigas = 'png';
bin = .02;
tmax = 1;
if ~isempty(varargin); assign(varargin{:}); end
day = idx(1,1); epoch = idx(1,2);
[~, fname,~] = fileparts(mfilename('fullpath'));
outdir = sprintf('%s/%s/%s/', pconf.andef{4},fname,animal,animal);

out.index = idx; out.data = []; out.dims = {};
% get the ripplekons events, which could be variably named
evid = find(contains(varargin(1:2:end), 'rippleskon'));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};
try 
    eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
catch; fprintf('no events detected for day%d ep%d\n', day,epoch); return; end
eventTime = eventTime(~isExcluded(eventTime(:,1),excludeIntervals),:);
swrStart = eventTime(:,1);

%% get licks
lickDIOIdx = task{day}{epoch}.inputdio;
isinput = cellfun(@(x) isequal(x.input,1), DIO{day}{epoch}, 'un', 1);
dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')), DIO{day}{epoch}, 'un', 1);
licks = sortrows(cell2mat(cellfun(@(x) [x.times repmat(str2double(regexp(x.original_id,'\d*','Match')), ...
    length(x.times),1)], DIO{day}{epoch}(ismember(dioID(isinput), lickDIOIdx)), 'un', 0)'),1);

%% plot
Pp=load_plotting_params({'defaults','dfa_lickswrcorr'});
if saveplots && ~displayplots; close all;
    ifig = figure('Visible','off','units','normalized','position', Pp.position, ...
        'color','white'); else
    ifig = figure('units','normalized','position',Pp.position, 'color','white'); end

% sf1 = subaxis(1,1,1,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
%     'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
%     'MarginBottom', Pp.MgBm); set(gca, 'Tag', 'ripkons');

realxcs = spikexcorr(swrStart,licks(:,1), bin, tmax);
bar(realxcs.time, realxcs.c1vsc2, 'k'); axis tight; hold on;
% time shift shuffle
x = [];
r = randi([-200 200],length(swrStart),1000)/1e3;
for i = 1:1000
    swrStartShift = sort(swrStart+r(:,i));
    xcs = spikexcorr(swrStartShift,licks(:,1), bin, tmax);
    x(i,:) = xcs.c1vsc2;
end
xstd = std(x);
xmean = mean(x);
plot(xcs.time, xmean, 'color', 'r', 'linewidth', 1);
fill([xcs.time'; flipud(xcs.time')],[xmean'-xstd';flipud(xmean'+xstd')],'b',...
    'linestyle','none', 'facealpha', .2);
ylim([min(realxcs.c1vsc2)-10 max(realxcs.c1vsc2)+10])
hold off
ylabel('xcorr')
xlabel('lag s')

% ---- super axis -----
sprtit = sprintf('%s %d %d lickVSswr bin20ms shuff200msSTD', animal,day,epoch);
sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center', ...
    'FontSize',12);
h = get(gcf,'Children');
set(gcf,'Children',flip(h)); % put super axis at bottom of axis stack. allows for zoom
% ---- pause, save figs ----
if displayplots; pause; end
if saveplots; save_figure(outdir, sprtit, 'savefigas', savefigas); end
end

