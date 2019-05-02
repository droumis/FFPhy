

%{

1. make swr iterator within analysis function
2. make resampling/KDE/binning operator for swr data chunks in analysis
function
3. make increment process creator to convert spiketimes to gridded logical
using a specified time vec as bin edges

4. combine spikes from same tet into mu spikes list

%FIRST.. just get xcorr measures for entire epochs (riptimes)

%}
%%

% ind = %day ep tetA clustA tetB clustB
% excludetimes = %include rips
% spikes = %spikes data for day (or just the two spiketrains)
% linpos = []; % won't do anything with default flags.. but is a required arg
% 
% calcxcorrmeasures(ind, excludetimes, spikes, linpos)
%%

%%
animals = {'JZ1'};
epochfilter = '(isequal($environment, ''wtrack''))';
tetfilter = '(~isequal($area, ''ca1''))';
% cellpairfilter = {'allcomb', '(~isequal($area, ''ca1'')', ...
%     '(~isequal($area, ''ca1''))'}; %need to make mua a field itself in cellinfo to make it easier to filter for
cellpairfilter = {'allcomb', ...
    '(isequal($area, ''mec'') && $numspikes>100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''),$tags,''un'',0)))))',...
    '(isequal($area, ''mec'') && $numspikes>100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''),$tags,''un'',0)))))'};

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 2;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
timefilter{1} = {'getconstimes', '($cons == 1)', ...
    'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
    'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
    'minvelocity', minvelocity,'maxvelocity',maxvelocity};

f = createfilter('animal', animals, 'days', [1:4], 'epochs', epochfilter, ...
    'cellpairs', cellpairfilter, 'iterator', 'singlecellanal', ...
    'excludetimefilter', timefilter);

f = setfilterfunction(f, 'calcxcorrmeasures', {'spikes', 'linpos'});
%%
fout = runfilter(f);

%%

fout.output{1}(1)

arrayfun(@(x) x.index(1:2) == [1 2], fout.output{1}, 'un', 0)
d = 1;
e = 2;


matInds = cell2mat({fout.output{1}.index}');
dayep_outinds = find(ismember(matInds(:,1:2), [d e], 'rows'));
coz = [fout.output{1}(dayep_outinds).coactivez]';
inds = cell2mat({fout.output{1}(dayep_outinds).index}');

allcoz = [fout.output{1}(:).coactivez]';
allinds = cell2mat({fout.output{1}(:).index}');

% accumarray([x(:),y(:)],values(:))
% heatmap

%alternative method of getting indices into fout for this d e
% fun = @(x) ismember(fout(1).output{1}(x).index(1:2), [d e], 'rows');
% inds = find(arrayfun(fun, 1:numel(fout(1).output{1}), 'un',1));

for i = size(fout.epochs{1},2)
    d = fout.epochs{1}(i,1);
    e = fout.epochs{1}(i,2);
    % for this d,e get all results from the t,c pairs
    
    data_inds_in_dayep = fout.data{1}{1};
    fout.output
end
%% use all spikes as MU

animals = {'D13','JZ1'};
epochfilter = '(isequal($environment, ''wtrack''))';
tetpairfilter = {'(isequal($area, ''mec''))', '(isequal($area, ''mec''))'};

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 2;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
timefilter{1} = {'getconstimes', '($cons == 1)', ...
    'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
    'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
    'minvelocity', minvelocity,'maxvelocity',maxvelocity};

% just use singletetrodeanal for the pair bc it will receive all the spikes
% for the day and the data indices for both tetrodes of the pair
f = createfilter('animal', animals, 'epochs', epochfilter, ...
    'tetrodepairs', tetpairfilter, 'iterator', 'singletetrodeanal', ...
    'excludetimefilter', timefilter);
% multitetrodeanal is not a tet pair iter.. it loads ALL ntrodes specified
% in the f.eegdata for each epoch.. so it would be used for instead to
% detect ripples where you have a handful of tets that you want to use at
% once, iterating over epochs
% tet pair filter is getting put into eegdata bc i'm using the
% eegtetrodepairs flag.. but this can't be used with singletetrodeanal, 
% which looks in .data.. 
f = setfilterfunction(f, 'mua_calcxcorrmeasures', {'spikes', 'linpos'});

fout = runfilter(f);

%%
animal = 'D13';
andef = animaldef(animal);
bs = load([andef{2}, animal, 'BehaveState.mat']);

%%
clf
allcoz = [fout(1).output{1}(:).ec]';
allinds = cell2mat({fout(1).output{1}(:).index}');
allinds = allinds(allinds(:,1)<8, :);
allcoz = allcoz(allinds(:,1)<8, :);

subplot(2,2,1)
% for d = unique(allinds(:,1))'
de = unique(allinds(:,1:2), 'rows');
pairs = unique(allinds(:,3:4), 'rows');
for tp = 1:length(pairs(:,1))
%         i = ismember(allinds(:,[1 2]), [d e], 'rows');
%         accray = accumarray([allinds(i,3),allinds(i,4)],allcoz(i));
    v = ismember(allinds(:,3:4), pairs(tp,:), 'rows');
    plot(allcoz(v))
    hold on
end
xlim([1 14])
set(gca,'xticklabel',[])
hold off
ylabel('excess corr')
xlabel('epoch (2/day)')

p = [];
dp = [];
for i = 1:size(de,1);
    bdei = ismember(bs.BehaveState.statespace.allbound(:,5:6), de(i,:), 'rows');
    b = bs.BehaveState.statespace.allbound(bdei,1);
    p = [p mean(b)];
    dp = [dp mean(diff(b))];
end
subplot(2,2,3)
[ax,h1,h2] = plotyy(1:length(p),p,1:length(dp),dp);
xlim([1 14])
xlabel('epoch (2/day)')
ylabel(ax(1),'prob correct')
ylabel(ax(2),'change')

%% compute correlation/regression 
% set(0,'defaultfigurecolor',[1 1 1])
subplot(2,2,2)
plot(allcoz(v),dp,'.', 'MarkerSize', 10,'Color', [.8 .3 0], 'MarkerEdgeColor', [.8 .3 0])
corrcoef(allcoz(v),dp,'rows','pairwise');
lsline
ylabel('change')
xlabel('excess corr')

subplot(2,2,4)
plot(allcoz(v),p,'.', 'MarkerSize', 10,'Color', [.3 .7 .8], 'MarkerEdgeColor', [.3 .7 .8])
lsline
ylabel('prob correct')
xlabel('excess corr')

%%
sr2p = [];
sr2v = [];
for tp = 1:length(pairs(:,1))
    v = ismember(allinds(:,3:4), pairs(tp,:), 'rows');
    [r2v,r2p] = corrcoef(allcoz(v),p,'rows','pairwise');
    sr2p = [sr2p r2p(1,2)];
    sr2v = [sr2v r2v(1,2)];
end

[Bp,Ip] = sort(sr2p, 'ascend');
psortedpairs = pairs(Ip,:);

subaxis(2,2,1)
for tp = 1:length(psortedpairs(:,1))
    v = ismember(allinds(:,3:4), psortedpairs(tp,:), 'rows');
    plot(allcoz(v))
    hold on
    pause
end
xlim([1 14])
set(gca,'xticklabel',[])
hold off
ylabel('excess corr')
xlabel('epoch (2/day)')

%%
clf
s = psortedpairs(:,1)';
t = psortedpairs(:,2)';
weights = B;
G = graph(s,t,weights);
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
plot(G,'LineWidth',LWidths)

%%

for tp = 1:length(psortedpairs(:,1))
    v = ismember(allinds(:,3:4), psortedpairs(tp,:), 'rows');
    subplot(2,2,1)
    plot(allcoz(v),dp,'.', 'MarkerSize', 10,'Color', [.8 .3 0], 'MarkerEdgeColor', [.8 .3 0])
    lsline
    ylabel('change')
    xlabel('excess corr')
    title(Bp(tp))

    subplot(2,2,3)
    plot(allcoz(v),p,'.', 'MarkerSize', 10,'Color', [.3 .7 .8], 'MarkerEdgeColor', [.3 .7 .8])
    lsline
    ylabel('prob correct')
    xlabel('excess corr')
    pause
end















