

%{
- calculate mu calcxcorrmeasures
- plot rip triggered spike rasters (with zscored) per epoch (or day?),
 along with tet pair corr and performance

%}

Fp = load_filter_params('mua_calcxcorrmeasures');
Fp.animals = {'D10', 'D13', 'JZ1', 'JZ3', 'JZ4'};
Fp.days = [1:7];

runFilterFramework = 0;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
load_other_data = 0;
plotstuff = 1;

savefigs = 0;
pausefigs = 1;

%% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp);
%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework == 1
    F = createfilter('animal', Fp.animals, 'days', Fp.days, 'epochs', ...
        Fp.epochfilter, 'tetrodepairs', Fp.tetpairfilter, 'iterator', ...
        Fp.iterator, 'excludetimefilter', Fp.timefilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes);
    for a = 1:length(F) % save filter detailes along with results
        F(a).datafilter_params = Fp;
    end
    F = runfilter(F);
end
%% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput == 1;
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave)
end
%% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput == 1;
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave);
end
%% Load and process other data
if load_other_data
    me = animaldef('demetris');
    
    % RIP TRIG load, get animals
    filename = 'D10-D13-JZ1-JZ3-JZ4_dfa_riptrigspikings_wtrack.mat';
    filfunc = 'dfa_riptrigspiking';
    riptrig_result = load(sprintf('%s/%s/%s', me{2}, filfunc, filename));
    riptrig_result = riptrig_result.F;
    riptrig_animals = {riptrig_result.animal};
    b = vertcat(riptrig_animals{:});
    riptrig_animals = b(:,1);
    
    %% RIP CORR get animals
    ripcorr_result = F;
    ripcorr_animals = {ripcorr_result.animal};
    b = vertcat(ripcorr_animals{:});
    ripcorr_animals = b(:,1);
    for ani = 1:length(ripcorr_animals)
        % get ripcorr data indices with excess correlation val into mat
        ripcorr_result_mat{ani} = cell2mat({ripcorr_result(ani).output{1}.index}');
        ripcorr_result_mat{ani} = [ripcorr_result_mat{ani} ...
            [ripcorr_result(ani).output{1}(:).ec]'];
        
        % get rip trig data indices
        riptrig_inds_mat{ani} = cell2mat({riptrig_result(ani).output{1}.index}');
        
        % PERFORMANCE (behave state) per day/ep
        % load performance and gather animals, indices of data
        % probs need to remake/save bstate into resultdir instead of ffdir
        filename = sprintf('%sBehaveState.mat',ripcorr_animals{ani});
        andef = animaldef(ripcorr_animals{ani});
        bstate = load(sprintf('%s/%s', andef{2}, filename));
        stspace{ani} = bstate.BehaveState.statespace.allbound;
        stspace_times{ani} = bstate.BehaveState.statespace.allepsMat(:,3);
        
        %% compute corr coef between rip corr and performance/change
        % for each nt pair, get perf mean and diff, out [day ep results]
        % the reason i'm doing this
        % here instead of in the pair figure loop is so i can compute the
        % corrcoefs and therefore have a template graph as background for
        % each pair edge in a plot.. like have all non active edges alpha
        % low, but visible so that the edge width of the active pair is put
        % into some context.
        de_behave{ani} = [];
        for ires = 1:length(ripcorr_result_mat{ani}(:,1));
            day = ripcorr_result_mat{ani}(ires,1);
            ep = ripcorr_result_mat{ani}(ires,2);
            ssinds = ismember(stspace{ani}(:,[5 6]), [day ep], 'rows');
            de_mean = mean(stspace{ani}(ssinds,1));
            de_meandiff = mean(abs(diff(stspace{ani}(ssinds,1))));
            de_behave{ani} = [de_behave{ani}; day ep de_mean de_meandiff];
        end
        ripcorr_result_mat{ani} = [ripcorr_result_mat{ani} de_behave{ani}];
        
        pairs = unique(ripcorr_result_mat{ani}(:,[3 4]), 'rows');
        ripcorrVperform{ani} = [];
        for ipair = 1:size(pairs,1)
            pinds = ismember(ripcorr_result_mat{ani}(:,[3 4]), pairs(ipair,:), 'rows');
            [bv,bp] = corrcoef(ripcorr_result_mat{ani}(pinds, 5),...
                ripcorr_result_mat{ani}(pinds, 8), ...
                'rows', 'pairwise');
            [dbv,dbp] = corrcoef(ripcorr_result_mat{ani}(pinds, 5),...
                ripcorr_result_mat{ani}(pinds, 9), ...
                'rows', 'pairwise');
            ripcorrVperform{ani} = [ripcorrVperform{ani}; pairs(ipair,:) bv(1,2) bp(1,2) dbv(1,2) ...
                dbp(1,2)];
        end
        
    end
end
%% ---------------- Plot ------------------------------------------------------
if plotstuff
    rows = 7;
    Pp = load_plotting_params('ripcorrVperform_wriptrig');
    %% Figure loop = for each nt pair in corr, plot corr/trig/perf per day/ep
    for ani = 1:length(ripcorr_animals)
                %% make template graph from corrcoef result
                %here
        pairs = unique(ripcorr_result_mat{ani}(:,[3 4]), 'rows');
        for ipair = 1:length(pairs(:,1));
            if savefigs && ~pausefigs;
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            %% plot rip trig
            nta = pairs(ipair,1);
            ntb = pairs(ipair,2);
            % get rip trig results for each ntrode, across all day/eps
            nta_pinds = find(ismember(riptrig_inds_mat{ani}(:,3), nta, 'rows'));
            ntb_pinds = find(ismember(riptrig_inds_mat{ani}(:,3), ntb, 'rows'));
            nta_dayeps = riptrig_inds_mat{ani}(nta_pinds,1:2);
            ntb_dayeps = riptrig_inds_mat{ani}(ntb_pinds,1:2);
            % day eps for both nt's should be the same if i make the rip trig
            % spiking and include the mu clusters..
            unq_dayeps = unique(nta_dayeps, 'rows');
            num_unq_dayeps = size(unq_dayeps,1)-5;
            for dei = 1:num_unq_dayeps
                day = unq_dayeps(dei,1);
                ep = unq_dayeps(dei,2);
                ntAriptrig_resultInds = find(ismember(riptrig_inds_mat{ani}(:, ...
                    [1 2 3]), [day ep nta], 'rows'));
                ntBriptrig_resultInds = find(ismember(riptrig_inds_mat{ani}(:, ...
                    [1 2 3]), [day ep ntb], 'rows'));
                % plot ntA riptrig day/eps and overlay clusters
                
                for i = 1:length(ntAriptrig_resultInds)
                    % plot rasters
                    subaxis(rows,num_unq_dayeps,dei)
                    cmap = colormap(lines);
                    di = riptrig_result(ani).output{1}(ntAriptrig_resultInds(i));
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
                    set(gca, 'YTickLabel',[])
%                     ylabel('rip num')
                    title(sprintf('%s.%d.%d.%d',riptrig_animals{ani},di.index(1:3)))
                    
                    %plot histogram (MAKE ZSCORE INSTEAD)
                    subaxis(rows,num_unq_dayeps,dei+num_unq_dayeps)
                    histogram(x1,60, 'FaceColor', cmap(i,:),...
                        'DisplayStyle', 'stairs');
                    hold on
                    h = histogram(x1,60, 'FaceColor', cmap(i,:), 'FaceAlpha', .2,...
                        'DisplayStyle', 'bar', 'EdgeAlpha', .1);
                    line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
                    set(gca, 'XTickLabel',[])
                    set(gca, 'YTickLabel',[])
                    axis tight
                    %     xlabel('time (s) from rip start')
%                     ylabel('spike count')
                    
                    %                     %% Create smoothed PSTH of mean firing rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % smoothedpsth = smoothvect(sum(ienvipsth,1)./(Fp.binsize*idata.noevents),Pp.kernel);
                    % ah = area(idata.time,zscore(smoothedpsth),'facecolor',...
                    %     [.2 .2 .2]);
                    %
                    % set(ah, 'FaceAlpha', Pp.areaAlpha)
                    % set(gca,'FontSize',8,'FontWeight','bold')
                end
                % plot ntB riptrig day/eps and overlay clusters
                for i = 1:length(ntBriptrig_resultInds)
                    cmap = colormap(lines);
                    di = riptrig_result(ani).output{1}(ntBriptrig_resultInds(i));
                    subaxis(rows,num_unq_dayeps,dei+2*num_unq_dayeps)
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
                    set(gca, 'YTickLabel',[])
%                     ylabel('rip num')
                    title(sprintf('%s-%d-%d-%d',riptrig_animals{ani},di.index(1:3)))
                    
                    subaxis(rows,num_unq_dayeps,dei+3*num_unq_dayeps)
                    histogram(x1,60, 'FaceColor', cmap(i,:),...
                        'DisplayStyle', 'stairs');
                    hold on
                    h = histogram(x1,60, 'FaceColor', cmap(i,:), 'FaceAlpha', .2,...
                        'DisplayStyle', 'bar', 'EdgeAlpha', .1);
                    line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
%                     set(gca, 'XTickLabel',round(di.time, 3))
                    set(gca, 'XTickLabel',[])
                    set(gca, 'YTickLabel',[])
                    axis tight
%                     xlabel('time (s) from rip start')
%                     ylabel('spike count')
                    hold on
                end
            end

            %% plot rip corr over day eps
            subaxis(rows,num_unq_dayeps,num_unq_dayeps*4+1:num_unq_dayeps*5)
            pinds = ismember(ripcorr_result_mat{ani}(:,[3 4]), pairs(ipair,:), 'rows');
            plot(ripcorr_result_mat{ani}(pinds,5));
            %% plot performance/change over day/eps    
            subaxis(rows,num_unq_dayeps,num_unq_dayeps*5+1:num_unq_dayeps*6)
            plot(ripcorr_result_mat{ani}(pinds,8))
            subaxis(rows,num_unq_dayeps,num_unq_dayeps*6+1:num_unq_dayeps*7)
            plot(ripcorr_result_mat{ani}(pinds,9))

            % % sslower = bs_statespace_allbound(:,2);
            % % sshigher = bs_statespace_allbound(:,3);
            % % % day = bs_statespace_allbound(:,5);
            % % % epoch = bs_statespace_allbound(:,6);
            % % errfillAll = fill([bs_statespace_time; flipud(bs_statespace_time)], ...
            % %     [sslower; flipud(sshigher)],[0 0 1],'linestyle','none');
            % % set(errfillAll, 'FaceAlpha', .1)
            % % hold on
            % % plot(bs_statespace_time, ssmode, '.')
            % % ylabel('statespace prob correct')            
            
            %% TODO: add rip trig lfp, occnorm wtrack firing, autocorr per day/ep
            % also add graph with the current pair edge highlighted.
            % seperately (beforehand).. make the graph primitive (with a
            % histology layer?) of the corrcoeff for all the pair ripcorr and the
            % performance
            
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s - nt%d nt%d', riptrig_animals{ani},nta,ntb);
            iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
                'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center');
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(paths.figdirectory, paths.filenamesave, sprtit)
            end
            close all
        end
    end
end

