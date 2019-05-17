
% plot the riptrig rasters for all the ntrodes together, over days

animals = {'D12'};
filtfunction = 'dfa_riptrigspiking';
env = 'wtrack';
me = animaldef('demetris');

loaddata = 1;
combine_per_ntrode = 1;
plotfigs = 1;
pausefigs = 0;
savefigs = 1;
justSU = 0;
justMEC = 1;
justMU = 0;

if ~exist('cmap') % only run tis once bc it launchs a figure bc matlab is stupip
    cmap = colormap(lines);
end
%%
% Fp = load_filter_params(filtfunction);
paths = make_paths(filtfunction, env);
if loaddata
    data = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, ...
        animals);
end

%% takes the stack of results per an,day,ep,nt,cl and returns an,nt with cl
% combined and epochs vertstacked. mu gets combined into one mucluster

if combine_per_ntrode
    psth = combine_riptrigspiking_perntrode(data);
end

%% plot each ntrode raster, mu and su specified. days as columns, ntrodes as rows
if plotfigs
    Pp = load_plotting_params('all_nts_days_riptrigspiking');
    for ian = 1:numel(animals)
        animal = animals{ian};
        fprintf('%s\n',animal);
        %% ---- init fig----
        if savefigs && ~pausefigs;
            close all
            ifig = figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            ifig = figure('units','normalized','position',Pp.position);
        end
        set(gcf,'color','white')
        %%
        allnts = [psth{ian}.ntrode];
        ntrodes = allnts;
        if justMEC
            %hack.. i should be doing a lookup into tetinfo
            % also, note, nt1 D10 is ca1. and JZ2 is all kinds of different
            ntrodes = ntrodes(ntrodes<16);
        end
        days = psth{ian}(1).days;
        numntrodes = numel(ntrodes);
        numdays = numel(days);
        ha = tight_subplot(numntrodes, numdays,[.01 .005],[.05 .05],[.05 .05]);
        for int = 1:numntrodes
            nt = ntrodes(int);
            ntind = find(allnts==nt);
%             ha = tight_subplot(numntrodes, numdays, .1, .1, .1);
            for idy = 1:numdays
%                 subplot(numntrodes, numdays, idy+((int-1)*numdays))
                axes(ha(idy+((int-1)*numdays)));
                if ~isempty(psth{ian}(ntind).mucluster)
                    day_mu_psth = cell2mat(psth{ian}(ntind).mucluster{days(idy)});
                    [sx, sy] = find(day_mu_psth');
                    f1 = scatter(sx,sy,2, 'filled');
                    %                     set(gca,'YDir','reverse')
                    if justSU
                        f1.MarkerFaceAlpha = 0; %hide mu
                    else
%                         f1.MarkerFaceAlpha = 0.02; 
                        f1.MarkerFaceAlpha = 0.1; %for JZ1 (sparser firing)
                    end
                    
                    f1.MarkerFaceColor = [0 0 0];
                    axis tight
                    set(gca,'XTick',[]);
                    set(gca,'YTick',[]);
                    if idy == 1;
                        ylabel(sprintf('%d',nt))
                    end
                    if int == numntrodes
                        xlabel(sprintf('%d',days(idy)))
                        set(gca,'XTick', [1 size(day_mu_psth,2)-1]);
                        if idy == 1;
                            set(gca,'XTickLabel', round(psth{ian}(ntind).psth_time([1, end]),1));
                        end
                    end
                    % ripple line
                    centerline = round(size(day_mu_psth,2)/2);
                    lh = line([centerline centerline], [1 size(day_mu_psth,1)]);
                    lh.LineWidth = 1;
                    lh.Color = [.8 0 0 .8];
                    % epoch lines
                    le = line([1 size(day_mu_psth,2)], ...
                        repmat(psth{ian}(ntind).eplengths{days(idy)}(1),2,1), 'Color', [0 0 0]);
                end
                if justMU
                    continue
                end
                if ~isempty(psth{ian}(ntind).suclusters)
                    try
                        if ~isempty(psth{ian}(ntind).suclusters{days(idy)})
                            hold on;
                            for isu = 1:length(psth{ian}(ntind).suclusters{days(idy)})
                                day_su_psth = psth{ian}(ntind).suclusters{days(idy)}{isu};
                                [sx, sy] = find(day_su_psth');
                                f1 = scatter(sx,sy,3, 'filled');
                                f1.MarkerFaceAlpha = .2;
                                f1.MarkerFaceColor = cmap(isu,:);
                                axis tight
                                set(gca,'XTick',[]);
                                set(gca,'YTick',[]);
                            end
                        end
                    catch
                    end
                end
                
            end
        end
        
        %% ---- super title and colorbar----
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
%         if justSU
         sprtit = sprintf('%s %s riptrigmu', env, animal); % just change the stitle manually for diff conditions
%         else
%             sprtit = sprintf('%s %s riptrigmusu', env, animal);
%         end
        iStitle = text(.5, .99, {sprtit}, 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center','FontSize', 16);
        
        spylabel = text(.01, .5, sprintf('ntrode'), 'Parent', sprtitleax, 'Units', ...
            'normalized', 'Rotation', 90);
        set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center', 'FontSize', 16);
        
        spxlabel = text(.5, .01, sprintf('day'), 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center', 'FontSize', 16);
        %% ---- pause, save figs ----
        if pausefigs
            pause
        end
        if savefigs
            save_figure(paths.figdirectory, paths.filenamesave, sprtit)
            close all
        end
    end
end

%% MOVED TO FUNCTION /home/droumis/Src/Matlab/filterframework_dr/Functions/combine_riptrigspiking_perntrode.m
% if 0
%     for ian = 1:numel(animals)
%         animal = animals{ian};
%         fprintf('%s\n',animal);
%         % get multiunit keys
%         mu_keys = evaluatefilter(cellinfo{ian}, 'isequal($tags, {''mua''})');
%         su_keys = evaluatefilter(cellinfo{ian}, '~isequal($tags, {''mua''})');
%         idata = data{ian}.F.output{1};
%         data_keys = cell2mat({idata.index}');
%         ntrodes = unique(data_keys(:,3));
%         for int = 1:numel(ntrodes) % ntrode
%             ntrode = ntrodes(int);
%             nt_data_keys = data_keys(data_keys(:,3)==ntrodes(int),:);
%             days = unique(nt_data_keys(:,1));
%             psth{ian}(int).ntrode = ntrodes(int);
%             psth{ian}(int).days = days;
%             for idy = 1:numel(days) % day
%                 day = days(idy);
%                 per_ep_mu_stack = {};
%                 nt_day_keys = nt_data_keys(nt_data_keys(:,1)==days(idy,:),:);
%                 % instead.. start with the mu/su loop instead, and
%                 % combine/stack eps..
%                 % find and combine cluster results and across eps
%                 nt_mu_keys = nt_day_keys(find(ismember(nt_day_keys, ...
%                     mu_keys, 'rows')),:);
%                 nt_su_keys = nt_day_keys(find(ismember(nt_day_keys, ...
%                     su_keys, 'rows')),:);
%                 if ~isempty(nt_mu_keys)
%                     epochs = unique(nt_mu_keys(:,2));
%                     for iep = 1:numel(epochs) % epoch
%                         per_ep_mu_stack{iep,1} = 0;
%                         epoch = epochs(iep);
%                         nt_ep_mu_keys = nt_mu_keys(nt_mu_keys(:,2)==epoch,:);
%                         nt_ep_mu_inds = find(ismember(data_keys, nt_ep_mu_keys, ...
%                             'rows'));
%                         for k = 1:numel(nt_ep_mu_inds) % MU cluster
%                             try
%                                 per_ep_mu_stack{iep,1} = per_ep_mu_stack{iep,1} + ...
%                                     idata(nt_ep_mu_inds(k)).psth;
%                             catch
%                                 fprintf('could not add mu spikes from %d %d %d %d, skipping\n',nt_ep_mu_keys(k,:));
%                                 continue
%                             end
%                         end
%                     end
%                     psth{ian}(int).eplengths{days(idy)} = cellfun(@(x) ...
%                         length(x(:,1)), per_ep_mu_stack,'un',1);
%                     psth{ian}(int).psth_time = idata(nt_ep_mu_inds(k)).time;
%                 end
%                 psth{ian}(int).mucluster{days(idy)} = {cell2mat(per_ep_mu_stack)};
%                 if ~isempty(nt_su_keys)
%                     sus = unique(nt_su_keys(:,4));
%                     for isu = 1:numel(sus) % SU clusters
%                         sui_keys = nt_su_keys(ismember(nt_su_keys(:,4),sus(isu),'rows'),:);
%                         nt_ep_su_inds = find(ismember(data_keys, sui_keys, ...
%                             'rows'));
%                         psth{ian}(int).suclusters{days(idy)}{isu} = cell2mat({idata(nt_ep_su_inds).psth}');
%                         psth{ian}(int).sucluster_id{days(idy)}{isu} = sus(isu);
%                     end
%                 end
%             end
%         end
%     end
% end
