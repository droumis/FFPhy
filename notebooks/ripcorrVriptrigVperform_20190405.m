

% step 1 make per tet rip trig firing across days
% step 2 combine all pairs
% step 3 load performance, calculate ripcorr V perform
% step 4 plot perform, change, riptrig

% load riptrig data, cellinfo
me = animaldef('demetris');
animals = {'D10', 'D13', 'JZ1', 'JZ3', 'JZ4'};
env = 'wtrack';
loaddata = 0;
plotfigs = 1;
savefigs = 1;
pausefigs = 0;

% load data for each animal and get info about indices and which are mu

if loaddata
    for an =1:length(animals);
        f = 'dfa_riptrigspiking';
        tmp = load(sprintf('%s/%s/%s_%s_%s.mat',me{2}, f, f, env, animals{an}));
        andata{an} = tmp.F.output{1};
        % get [day ep tet cl] indices for the data output struct array
        data_indices{an} = cell2mat({andata{an}.index}');
        % get cell info
        andef{an} = animaldef(animals{an});
        cellinfo{an} = loaddatastruct(andef{an}{2}, andef{an}{3}, 'cellinfo');
        % get multiunit inds
        mu_ids{an} = evaluatefilter(cellinfo{an}, 'isequal($tags, {''mua''})');
        % get mapping to each unique tetrode
        [uniqtet{an},~, data2uniqtet{an}] = unique(data_indices{an}(:,3));
    end
    color = colormap(lines);
end 

% add titles per subplot,
% ripple line
% subplots tight
% histogram + zscore in second row
% add super title
% add figure saving
% make into functions so i can call by pairs to plot

paths = make_paths('riptrigspiking', 'wtrack');

if plotfigs
    Pp = load_plotting_params('riptrigspiking');
for an = 1:length(animals);
for i = 1:length(uniqtet{an})
    if savefigs && ~pausefigs;
        ifig = figure('Visible','off','units','normalized','position', ...
            Pp.position);
    else
        ifig = figure('units','normalized','position',Pp.position);
    end
    set(gcf,'color','white')
    tet_datainds = find(data2uniqtet{an} == i);
    tet_ids = data_indices{an}(tet_datainds,:);
    days = unique(tet_ids(:,1));
    for d = 1:length(days);
        % find this day for tet inds into data
        idt_datainds = find(tet_ids(:,1)==days(d));
        detc_ids = tet_ids(idt_datainds,:);

        % find and combine multiunit cluster results and acros eps
        detc_mu_ids = detc_ids(find(ismember(detc_ids,mu_ids{an},'rows')),:);
        if ~isempty(detc_mu_ids)
            % for each ep (combine mu clusters within ep)
            eps = unique(detc_mu_ids(:,2));
            mu_eps_psth = [];
            for e = 1:length(eps)
                mu_inep_ids = detc_mu_ids(detc_mu_ids(:,2)==eps(e),:);
                mu_inds = find(ismember(data_indices{an}, mu_inep_ids,'rows'));
                mu_psth = zeros(size(andata{an}(mu_inds(1)).psth));
                for m = 1:size(mu_inds,1);
                    mu_psth = mu_psth+andata{an}(mu_inds(m)).psth;
                end
                % stack eps vertically
                mu_eps_psth = [mu_eps_psth;mu_psth];
            end
            subaxis(2,length(days),d)
            hold on
            [x1, y1] = find(mu_eps_psth);
            f1 = scatter(x1,y1,2, 'filled');
            f1.MarkerFaceAlpha = 0.08;
            f1.MarkerFaceColor = [0 0 0];
            axis tight
            
            subaxis(2,length(days),d+length(days))
            hold on
            nbins=20;
%             histogram(x1,nbins, 'Color', [.8 .8 .8], 'DisplayStyle', 'stairs', 'Normalization','probability');
            
            histogram(x1,nbins, 'FaceColor', [.8 .8 .8], 'FaceAlpha', .1,...
                'DisplayStyle', 'bar', 'EdgeAlpha', .3, 'Normalization','probability');
        end

        detc_su_ids = detc_ids(find(~ismember(detc_ids,mu_ids{an},'rows')),:);
        if ~isempty(detc_su_ids)
            % for each unique SU on this tet, vert stack across eps
            [uniq_su_ids, ~, full2uniq_su] = unique(detc_su_ids(:,4));
            for s = 1:size(uniq_su_ids,1)
                su_ids_xeps = detc_su_ids(detc_su_ids(:,4) == uniq_su_ids(s),:);
                eps = unique(su_ids_xeps(:,2));
                su_eps_psth = [];
                for e = 1:length(eps)
                    su_iep_ind = find(ismember(data_indices{an}, su_ids_xeps(e,:), 'rows'));
                    su_psth = andata{an}(su_iep_ind).psth;
                    if sum(sum(su_psth)) < 20
                        continue
                    end
                    % stack eps vertically
                    su_eps_psth = [su_eps_psth; su_psth];
                end

                subaxis(2,length(days),d)
                [x1, y1] = find(su_eps_psth);
                f1 = scatter(x1,y1,3, 'filled');
                f1.MarkerFaceAlpha = .8;
                f1.MarkerFaceColor = color(s,:);

                %lower aggregate plots
                subaxis(2,length(days),d+length(days))
                % TODO normalize histograms 0-1 first
%                 histogram(x1,nbins, 'Color', color(s,:), 'DisplayStyle', 'stairs', 'Normalization','probability');
                hold on
                histogram(x1,nbins, 'FaceColor', color(s,:), 'FaceAlpha', .5,...
                    'DisplayStyle', 'bar', 'EdgeAlpha', .3, 'Normalization','probability');

            end
        end
        
        subaxis(2,length(days),d)
        ax = gca;
        lx = round(ax.XLim(end)/2);
        line([lx lx],ax.YLim,'Color',[0 0 0], 'LineStyle', '--')
        set(gca, 'XTickLabel',[])
        set(gca, 'YTickLabel',[])
        title(sprintf('d%d',days(d)));

        pos = ax.Position; % [x y width height]
        pos(1) = pos(1) - 0.02;
        pos(3) = 0.11;
        set(gca, 'Position', pos)
        axis tight
        
        subaxis(2,length(days),d+length(days))
        ax = gca;
        line([lx lx],ax.YLim,'Color',[0 0 0], 'LineStyle', '--')
        set(gca, 'XTickLabel',[])
        set(gca, 'YTickLabel',[])

        pos = ax.Position; % [x y width height]
        pos(1) = pos(1) - 0.02;
        pos(3) = 0.11;
        set(gca, 'Position', pos)
        axis tight
        
    end %days
    %% ---- super title and colorbar----
    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('riptrig MU-SU %s %s - nt%d', env, animals{an}, tet_ids(1,3));
    iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
        'normalized'); 
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center');
    Sylabel1 = text(.06, .7, 'swr #', 'Parent', sprtitleax, 'Units', ...
        'normalized'); 
    set(Sylabel1,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center', 'Rotation',90);
    Sylabel2 = text(.06, .3, 'spike probability', 'Parent', sprtitleax, 'Units', ...
        'normalized', 'Rotation',90); 
    set(Sylabel2,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
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