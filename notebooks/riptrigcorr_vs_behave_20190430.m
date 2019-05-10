


Fp = load_filter_params('mua_calcxcorrmeasures');
Fp.animals = {'D10', 'D13', 'JZ1', 'JZ3', 'JZ4'};
Fp.days = 1:7;

%%
F = createfilter('animal', Fp.animals, 'days', Fp.days, 'epochs', ...
    Fp.epochfilter, 'tetrodepairs', Fp.tetpairfilter, 'iterator', ...
    Fp.iterator, 'excludetimefilter', Fp.timefilter);
F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes);
for i = 1:length(F)
    F(i).datafilter_params = Fp;
end
%%
fout = runfilter(F);
ripcorr_result = fout;
%%
paths = make_paths(Fp);
%%
save_filter_output(ripcorr_result, paths.filtOutputDirectory, paths.filenamesave)
%% get ripcorr result's animals and indices

ripcorranimals = {ripcorr_result.animal};
b = vertcat(ripcorranimals{:});
ripcorranimals = b(:,1);

for ani = 1:length(ripcorranimals)
    ripcorr_result_mat{ani} = cell2mat({ripcorr_result(ani).output{1}.index}');
    ripcorr_result_mat{ani} = [ripcorr_result_mat{ani} [ripcorr_result(ani).output{1}(:).ec]'];
end

%% load riptrig spiking results and get animals, indices
filename = 'D10_dfa_riptrigspikings_openfield.mat';
deme = animaldef('demetris');
filfunc = 'dfa_riptrigspiking';
riptrig_result = load(sprintf('%s/%s/%s', deme{2}, filfunc, filename));
riptrig_result = riptrig_result.F;
riptriganimals = {riptrig_result.animal};
b = vertcat(riptriganimals{:});
riptriganimals = b(:,1);

for ani = 1:length(riptriganimals)
    riptrig_inds_mat{ani} = cell2mat({riptrig_result(ani).output{1}.index}');
%     riptrig_result_mat{ani} = [riptrig_result_mat{ani} [riptrig_result(ani).output{1}(:).psth]'];
end

%% load behave state
for ani = 1:length(ripcorranimals)
bstate = load(sprintf(...
    '/data2/demetris/%s/filterframework/%sBehaveState.mat',ripcorranimals{ani},...
    ripcorranimals{ani}));
stspace{ani} = bstate.BehaveState.statespace.allbound;
stspace_times{ani} = bstate.BehaveState.statespace.allepsMat(:,3);
end
%% get behave state mean ep val and diff for each record in results and concat
for ani = 1:length(ripcorranimals)
    de_behave{ani} = [];
    for ipair = 1:length(ripcorr_result_mat{ani}(:,1));
        d = ripcorr_result_mat{ani}(ipair,1);
        e = ripcorr_result_mat{ani}(ipair,2);
        ssinds = ismember(stspace{ani}(:,[5 6]), [d e], 'rows');
        de_mean = mean(stspace{ani}(ssinds,1));
        de_meandiff = mean(abs(diff(stspace{ani}(ssinds,1))));
        de_behave{ani} = [de_behave{ani}; d e de_mean de_meandiff];
    end
    ripcorr_result_mat{ani} = [ripcorr_result_mat{ani} de_behave{ani}];
end

%% plot the pair rip corr across epochs
pair_ripcorr = [];
for ani = 3%1:length(animals)
    pair_ripcorr{ani} = [];
    pairs = unique(ripcorr_result_mat{ani}(:,[3 4]), 'rows');
    for ipair = 1:length(pairs(:,1));
        pinds = ismember(ripcorr_result_mat{ani}(:,[3 4]), pairs(ipair,:), 'rows');
        pair_ripcorr{ani} = [pair_ripcorr{ani}; ripcorr_result_mat{ani}(pinds,5)'];
    end
    subaxis(4,1,[1 2])
    imagesc(pair_ripcorr{ani})
    xlabel('epoch')
    ylabel('pair ripcorr')
    xlim([1 14])
    title(ripcorranimals{ani})
    set(gca,'YTickLabel',[]);
    
    subaxis(4,1,3)
    plot(1:length(find(pinds)),ripcorr_result_mat{ani}(find(pinds),8));
    xlim([1 14])
    ylabel('prob correct')
    
    subaxis(4,1,4)
    plot(2:length(find(pinds))+1,ripcorr_result_mat{ani}(find(pinds),9));
    xlim([1 14])
    ylabel('change')
    xlabel('epoch')
    set(gca,'YTickLabel',[]);

    return
    %% save fig
    
    
end
set(gcf,'color','w');

%% Plot Pair rip trig over days.. TODO: ADD the ripCORR overdays
% currently waiting on wtrack riptrigspiking to run

pairs = unique(ripcorr_result_mat{ani}(:,[3 4]), 'rows');
for ipair = 1:length(pairs(:,1));
    nta = pairs(ipair,1);
    ntb = pairs(ipair,2);
    
    nta_pinds = find(ismember(riptrig_inds_mat{ani}(:,3), nta, 'rows'));
    ntb_pinds = find(ismember(riptrig_inds_mat{ani}(:,3), ntb, 'rows'));
    nta_dayeps = riptrig_inds_mat{ani}(nta_pinds,1:2);
    ntb_dayeps = riptrig_inds_mat{ani}(ntb_pinds,1:2);
    % day eps for both nt's should be the same if i make the rip trig
    % spiking and include the mu clusters
    for dei = 1:size(nta_dayeps, 2)
        d = nta_dayeps(dei,1);
        e = nta_dayeps(dei,2);
        ntAriptrig_resultInd = find(ismember(riptrig_inds_mat{ani}(:,[1 2 3]), ...
            [d e nta], 'rows'));
        ntBriptrig_resultInd = find(ismember(riptrig_inds_mat{ani}(:,[1 2 3]), ...
            [d e ntb], 'rows'));
        for i = 1:length(ntAriptrig_resultInd)
            cmap = colormap(lines);
            di = riptrig_result(ani).output{1}(ntAriptrig_resultInd(i));
            figure(4)
            subaxis(4,size(nta_dayeps,2),dei)
            hold on
            [y1, x1] = find(di.psth);
            f1 = scatter(x1,y1,3, 'filled');
            f1.MarkerFaceAlpha = 0.4;
            f1.MarkerFaceColor = cmap(i,:);
            axis tight
            lx = round(size(di.psth,2)/2);
            ly = size(di.psth,1);
            line([lx lx],[1 ly],'Color',[0 0 0], 'LineStyle', '--')
            set(gca, 'XTickLabel',[])
            ylabel('rip num')
            title(sprintf('%s-%d-%d-%d',riptriganimals{ani},di.index(1:3)))

            subaxis(4,size(nta_dayeps,2),dei+size(nta_dayeps,2))
            histogram(x1,60, 'FaceColor', cmap(i,:),...
                'DisplayStyle', 'stairs');
            hold on 
            h = histogram(x1,60, 'FaceColor', cmap(i,:), 'FaceAlpha', .2,...
                'DisplayStyle', 'bar', 'EdgeAlpha', .1);
            line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
            set(gca, 'XTickLabel',[])
            axis tight
        %     xlabel('time (s) from rip start')
            ylabel('spike count')
            hold on
        end
        for i = 1:length(ntBriptrig_resultInd)
            cmap = colormap(lines);
            di = riptrig_result(ani).output{1}(ntBriptrig_resultInd(i));
            figure(3)
            subaxis(4,size(nta_dayeps,2),dei+2*size(nta_dayeps,2))
            hold on
            [y1, x1] = find(di.psth);
            f1 = scatter(x1,y1,3, 'filled');
            f1.MarkerFaceAlpha = 0.4;
            f1.MarkerFaceColor = cmap(i,:);
            axis tight
            lx = round(size(di.psth,2)/2);
            ly = size(di.psth,1);
            line([lx lx],[1 ly],'Color',[0 0 0], 'LineStyle', '--')
            set(gca, 'XTickLabel',[])
            ylabel('rip num')
            title(sprintf('%s-%d-%d-%d',riptriganimals{ani},di.index(1:3)))

            subaxis(4,size(nta_dayeps,2),dei+3*size(nta_dayeps,2))
            histogram(x1,60, 'FaceColor', cmap(i,:),...
                'DisplayStyle', 'stairs');
            hold on 
            h = histogram(x1,60, 'FaceColor', cmap(i,:), 'FaceAlpha', .2,...
                'DisplayStyle', 'bar', 'EdgeAlpha', .1);
            line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
            set(gca, 'XTickLabel',round(di.time, 3))
            axis tight
            xlabel('time (s) from rip start')
            ylabel('spike count')
            hold on 
        end
    end
    return
end

%%  plot the rip trig for each epoch (horiz) for each pair (vertstacked)
d = 1;
e = 6;
ntA = 6;
ntB = 11;
clf
ntAriptrig_resultInds = find(ismember(riptrig_inds_mat{ani}(:,[1 2 3]), ...
    [d e ntA], 'rows'));
ntBriptrig_resultInds = find(ismember(riptrig_inds_mat{ani}(:,[1 2 3]), ...
    [d e ntB], 'rows'));
for i = 1:length(ntAriptrig_resultInds)
    cmap = colormap(lines);
    di = riptrig_result(ani).output{1}(ntAriptrig_resultInds(i));
    figure(3)
    subaxis(4,1,1)
    hold on
    [y1, x1] = find(di.psth);
    f1 = scatter(x1,y1,3, 'filled');
    f1.MarkerFaceAlpha = 0.4;
    f1.MarkerFaceColor = cmap(i,:);
    axis tight
    lx = round(size(di.psth,2)/2);
    ly = size(di.psth,1);
    line([lx lx],[1 ly],'Color',[0 0 0], 'LineStyle', '--')
    set(gca, 'XTickLabel',[])
    ylabel('rip num')
    title(sprintf('%s-%d-%d-%d',riptriganimals{ani},di.index(1:3)))

    subaxis(4,1,2)
    histogram(x1,60, 'FaceColor', cmap(i,:),...
        'DisplayStyle', 'stairs');
    hold on 
    h = histogram(x1,60, 'FaceColor', cmap(i,:), 'FaceAlpha', .2,...
        'DisplayStyle', 'bar', 'EdgeAlpha', .1);
    line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
    set(gca, 'XTickLabel',[])
    axis tight
%     xlabel('time (s) from rip start')
    ylabel('spike count')
    hold on
end
for i = 1:length(ntBriptrig_resultInds)
    cmap = colormap(lines);
    di = riptrig_result(ani).output{1}(ntBriptrig_resultInds(i));
    figure(3)
    subaxis(4,1,3)
    hold on
    [y1, x1] = find(di.psth);
    f1 = scatter(x1,y1,3, 'filled');
    f1.MarkerFaceAlpha = 0.4;
    f1.MarkerFaceColor = cmap(i,:);
    axis tight
    lx = round(size(di.psth,2)/2);
    ly = size(di.psth,1);
    line([lx lx],[1 ly],'Color',[0 0 0], 'LineStyle', '--')
    set(gca, 'XTickLabel',[])
    ylabel('rip num')
    title(sprintf('%s-%d-%d-%d',riptriganimals{ani},di.index(1:3)))

    subaxis(4,1,4)
    histogram(x1,60, 'FaceColor', cmap(i,:),...
        'DisplayStyle', 'stairs');
    hold on 
    h = histogram(x1,60, 'FaceColor', cmap(i,:), 'FaceAlpha', .2,...
        'DisplayStyle', 'bar', 'EdgeAlpha', .1);
    line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
    set(gca, 'XTickLabel',round(di.time, 3))
    axis tight
    xlabel('time (s) from rip start')
    ylabel('spike count')
    hold on 
end

%% plot rip trig scratch
subaxis(6,1,[1 2])
[x1,y1] = find(a.F.output{1}(1).psth);
scatter(x1,y1, 'k.')
axis tight
SP = 134;
line([SP SP],[1 1000],'Color',[1 0 0])
set(gca, 'XTickLabel',[])
ylabel('rip num')
title(sprintf('%s-%d-%d-%d-%d',ripcorranimals{ani},a.F.output{1}(1).index))

subaxis(6,1,3)
histogram(x1,60, 'FaceColor', [0 0 0])
line([SP SP],[1 50],'Color',[1 0 0])
axis tight
xlabel('time (s) from rip start')
ylabel('spike count')

subaxis(6,1,[4 5])
[x1,y1] = find(a.F.output{1}(2).psth);
scatter(x1,y1, 'k.')
axis tight
SP = 134;
line([SP SP],[1 1000],'Color',[1 0 0])
set(gca, 'XTickLabel',[])
ylabel('rip num')
title(sprintf('%s-%d-%d-%d-%d',ripcorranimals{ani},a.F.output{1}(2).index))

subaxis(6,1,6)
histogram(x1,60, 'FaceColor', [0 0 0])
line([SP SP],[1 50],'Color',[1 0 0])
axis tight
xlabel('time (s) from rip start')
ylabel('spike count')

%% for each pair (3 4), compute cc of pair ripcorr (5) vs behave cols (8,9)
% put the cc into a new mat with the pair ID's 



for ani = 1:length(ripcorranimals);
    outbv{ani} = [];
    pairs = unique(ripcorr_result_mat{ani}(:,[3 4]), 'rows');
    for ipair = 1:length(pairs(:,1))
        pinds = ismember(ripcorr_result_mat{ani}(:,[3 4]),...
            pairs(ipair,:), 'rows');
        [bv,bp] = corrcoef(ripcorr_result_mat{ani}(pinds, 5),...
            ripcorr_result_mat{ani}(pinds, 8), ...
            'rows', 'pairwise');
        [dbv,dbp] = corrcoef(ripcorr_result_mat{ani}(pinds, 5),...
            ripcorr_result_mat{ani}(pinds, 9), ...
            'rows', 'pairwise');
        outbv{ani} = [outbv{ani}; pairs(ipair,:) bv(1,2) bp(1,2) dbv(1,2) ...
            dbp(1,2)];
    end
end

[Bp,Ip] = sort(outbv{3}(:,6), 'ascend');
outbv{3}(Ip,:)







