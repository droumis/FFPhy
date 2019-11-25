
% get occ norm firing rate of all single units

run_ff = 0;
savedata = run_ff;

load_result = 0;
% stack_lfp = 1;
% load_events = 1;
% make_filter_vecs = 1;

plotfigs = 1;
plot_wtrack_SFRmaps = 0;
plot_openfield_SFRmaps = 1;

savefigs = 0;
pausefigs = 1;

% use_filters = {'firstwell', 'noise'};
animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ4'}; %, 'D13'}; % next run ff jz1, jz3, jz4
add_params = {'wtrack'};

Fp.animals = animals;
Fp.add_params = add_params;
Fp.filtfunction = 'dfa_occNormFiring';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
filetail = '';

%% run filter/func
if run_ff == 1
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'cells',Fp.cellfilter,...
        'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F);
    for d = 1:length(F)
        F(d).datafilter_params = Fp;
    end
    %% stack result before saving
    for ian = 1:length(F)
        stacked = {};
        for de = 1:length(F(ian).output{1})
            stacked{de,1} = [F(ian).output{1}{de}{:}];
        end
        F(ian).outputstack = [stacked{:}];
    end
%     if run_Lin
%         lF = createfilter('animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
%         lF = setfilterfunction(lF, 'dfa_filtercalclinfields', {'spikes', 'linpos', 'cellinfo'});
%         lF = runfilter(lF);
%     end
    
end
%% save data
if savedata == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s%s', Fp.epochEnvironment, filetail))
end
%% load result
if load_result
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', Fp.epochEnvironment, filetail));
    openF = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', 'openfield', filetail));
end
%% plot
if plotfigs
    if plot_wtrack_SFRmaps
        Pp = load_plotting_params({'defaults', 'frmaps'});
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp(arrayfun(@(x) F(x).animal{3}, 1:length(F), 'un', 0), animal));
            data_keys = cell2mat({F(ian).outputstack.index}');
            opendata_keys = cell2mat({openF(ian).outputstack.index}');
            % cells are defined as per day, tet, cluster
            [unqclusts, ~, unqclust_idx] = unique(data_keys(:,[1 3 4]), 'rows');
            %% ---------- Get tet, task, cell tags --------------------------------
            andef = animaldef(animal);
            cellinfo = loaddatastruct(andef{2}, animal, 'cellinfo');
            tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
            taskinfo = loaddatastruct(andef{2}, animal, 'task');
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            pos = loaddatastruct(andef{2}, animal, 'pos');
            for cl = 1:length(unqclusts(:,1)) % plot per cluster
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                    ifig.InvertHardcopy = 'off';
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                day = unqclusts(cl,1);
                ntrode = unqclusts(cl,2);
                cluster = unqclusts(cl,3);
                cl_data_idx = find(ismember(data_keys(:,[1 3 4]), unqclusts(cl,:), 'rows'));
                
                numeps = length(cl_data_idx);
                fr = {};
                for e = 1:numeps
                    epoch = data_keys(cl_data_idx(e),2);
                    num_trajs = length(F(ian).outputstack(cl_data_idx(e)).smoothedspikerate);
                    for j = 1:num_trajs
                        % plot trajs (rows) x eps (cols)
                        %                         if numeps < 2
                        %                             rows = 2;
                        %                         else
                        %                             rows = numeps;
                        %                         end
                        sf = subaxis(2,num_trajs+2, j+(num_trajs+2).*(e-1), 'SpacingVert', Pp.SpVt, ...
                            'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                            Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                        
                        xt = F(ian).outputstack(cl_data_idx(e)).xticks{j};
                        yt = F(ian).outputstack(cl_data_idx(e)).yticks{j};
                        fr{end+1} = [];
                        try
                            fr{end}(xt,yt) = F(ian).outputstack(cl_data_idx(e)).smoothedspikerate{j};
                        catch
                            continue
                        end
                        fr{end} = fr{end}';
                        Zfr = reshape( zscore(fr{end}(:)), size(fr{end}) );
                        % plot background all pos
                        
                        plot(pos{day}{epoch}.data(:,6), pos{day}{epoch}.data(:,7), ...
                            'Color', [.2 .2 .2], 'LineWidth', 1); %[0 .1 .8], 
                        axis tight
                        xlimx = xlim;
                        ylimy = ylim;
                        hold on
                        B = imgaussfilt(Zfr, 2);
                        s = pcolor(B);
                        newmap = jet; %brighten(jet,-.2); %flipud(hot); %brighten(hot,.2));
                        colormap(gca,newmap)
                        mask1 = B > std2(B); % baseline empty val transparent
                        %caxis([0 alltrajPeakFR]) %max(max((I)))-std2(I)]) % saturate top standard dev
                        s.AlphaData = mask1; % make all zero bins transparent
                        s.EdgeColor = 'none';
                        s.FaceAlpha = 'interp';
                        s.FaceColor = 'interp';
                        xlim(xlimx)
                        ylim(ylimy)
                        set(sf, 'color', 'black')
                        hold off
                        xlabel('cm', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                            Pp.FontNm)
                        ylabel('cm','FontSize',Pp.FontS, 'FontWeight',Pp.FontW,'FontName', ...
                            Pp.FontNm)
                        if ~(j == 1 && e ~=1)
                            xlabel('')
                            ylabel('')
                            set(gca, 'ytick', []);
                            set(gca, 'xtick', []);
                        end
                        if e == 1
                            title(sprintf('%s', Pp.tits{j}), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                                'FontName', Pp.FontNm)
                        end
                    end
                    sf = subaxis(2,num_trajs+2, [5 6 13 14], 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    dims = cell2mat(cellfun(@size, fr,'un',0)');
                    maxx = max(dims(:,1));
                    maxy = max(dims(:,2));
                    frmapmean = zeros(maxx, maxy, num_trajs*numeps);
                    for i = 1:length(fr)
                        [xs, ys] = size(fr{i});
                        % insert into z slize of a
                        frmapmean(1:xs,1:ys,i) = fr{i};
                    end
                    frmapmean = squeeze(mean(frmapmean, 3));
                    plot(pos{day}{epoch}.data(:,6), pos{day}{epoch}.data(:,7), ...
                        'Color', [.2 .2 .2], 'LineWidth', 1);
                    axis tight
                    xlimx = xlim;
                    ylimy = ylim;
                    hold on
                    
                    frmapmeanZ = reshape( zscore(frmapmean(:)), size(frmapmean) );
                    Bm = imgaussfilt(frmapmeanZ, 2);
                    s = pcolor(Bm);
                    newmap = jet; %brighten(jet,-.2); %flipud(hot); %brighten(hot,.2));
                    colormap(gca,newmap)
                    mask1 = Bm > std2(Bm); % baseline empty val transparent
                    %caxis([0 alltrajPeakFR]) %max(max((I)))-std2(I)]) % saturate top standard dev
                    s.AlphaData = mask1; % make all zero bins transparent
                    s.EdgeColor = 'none';
                    s.FaceAlpha = 'interp';
                    s.FaceColor = 'interp';
                    xlim(xlimx)
                    ylim(ylimy)
                    set(sf, 'color', 'black')
                    hold off
                    xlabel('')
                    ylabel('')
                    set(gca, 'ytick', []);
                    set(gca, 'xtick', []);
                    title(sprintf('all eps trajs', Pp.tits{j}), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    clrbar = colorbar;
                    ylabel(clrbar, 'zscore firing rate');
                    
                    %% plot openfield
%                     sf = subaxis(2,num_trajs+2, [7 8 15 16], 'SpacingVert', Pp.SpVt, ...
%                         'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
%                         Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    % concat, mean, all open field epochs 
%                     xt = openF(ian).outputstack(cl_data_idx(e)).xticks{j};
%                     yt = openF(ian).outputstack(cl_data_idx(e)).yticks{j};
%                     fr{end+1} = [];
%                     fr{end}(xt,yt) = F(ian).outputstack(cl_data_idx(e)).smoothedspikerate{j};
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %d (%s) %d %d %s %s%s',Fp.paths.filenamesave(5:end), animal, ...
                    day, num2str(data_keys(cl_data_idx(:),2)'), ntrode, cluster, Fp.epochEnvironment, filetail, datestr(now, 'yyyymmdd'));
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure(Fp.paths.figdirectory, Fp.paths.filenamesave, sprtit)
                    close all
                end
            end
        end
    end
    if plot_openfield_SFRmaps
        Pp = load_plotting_params({'defaults', 'frmaps'});
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp(arrayfun(@(x) openF(x).animal{3}, 1:length(openF), 'un', 0), animal));
            data_keys = cell2mat({openF(ian).outputstack.index}');
            % cells are defined as per day, tet, cluster
            [unqclusts, ~, unqclust_idx] = unique(data_keys(:,[1 3 4]), 'rows');
            %% ---------- Get tet, task, cell tags --------------------------------
            andef = animaldef(animal);
            cellinfo = loaddatastruct(andef{2}, animal, 'cellinfo');
            tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
            taskinfo = loaddatastruct(andef{2}, animal, 'task');
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            pos = loaddatastruct(andef{2}, animal, 'pos');
            for cl = 1:length(unqclusts(:,1)) % plot per cluster
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                    ifig.InvertHardcopy = 'off';
                else
                    ifig = figure(1);%'units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                day = unqclusts(cl,1);
                ntrode = unqclusts(cl,2);
                cluster = unqclusts(cl,3);
                cl_data_idx = find(ismember(data_keys(:,[1 3 4]), unqclusts(cl,:), 'rows'));
                idx = [];
                poseps = [];
                for i = 1:length(cl_data_idx)
                    xt{i} = openF(ian).outputstack(cl_data_idx(i)).xticks{1};
                    yt{i} = openF(ian).outputstack(cl_data_idx(i)).yticks{1};
                    idx(i,:) = openF(ian).outputstack(cl_data_idx(i)).index;
                    poseps = [poseps; pos{day}{idx(i,2)}.data(:,6), pos{day}{idx(i,2)}.data(:,7)];
                end         
                dims = [cell2mat(cellfun(@length , xt, 'un', 0)') cell2mat(cellfun(@length , yt, 'un', 0)')];
                maxx = max(dims(:,1));
                maxy = max(dims(:,2));
                minx = max(dims(:,1));
                miny = max(dims(:,2));
                d = zeros(maxx, maxy, length(cl_data_idx));
                for i = 1:length(cl_data_idx)
                    try
                        d(xt{i},yt{i},i) = openF(ian).outputstack(cl_data_idx(i)).smoothedspikerate{1};
                    catch
                        continue
                    end
                end
                d = squeeze(mean(d,3));

                sf = subaxis(2,2, 1, 'SpacingVert', Pp.SpVt, ...
                    'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                    Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                
                    plot(poseps(:,1),poseps(:,2), 'Color', [.2 .2 .2], 'LineWidth', 1);
                    axis tight
                    xlimx = xlim;
                    ylimy = ylim;
                    hold on
                    
                    Z = reshape( zscore(d(:)), size(d));
                    Bm = imgaussfilt(Z, 2);
                    s = pcolor(xt{i}, 1:maxy, Bm);
                    newmap = jet; %brighten(jet,-.2); %flipud(hot); %brighten(hot,.2));
                    colormap(gca,newmap)
                    mask1 = Bm > std2(Bm)/10; % baseline empty val transparent
                    %caxis([0 alltrajPeakFR]) %max(max((I)))-std2(I)]) % saturate top standard dev
                    s.AlphaData = mask1; % make all zero bins transparent
                    s.EdgeColor = 'none';
                    s.FaceAlpha = 'interp';
                    s.FaceColor = 'interp';
                    xlim(xlimx)
                    ylim(ylimy)
                    set(sf, 'color', 'black')
                    hold off
                    xlabel('')
                    ylabel('')
                    set(gca, 'ytick', []);
                    set(gca, 'xtick', []);
                    title(sprintf('all eps trajs', Pp.tits{j}), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    clrbar = colorbar;
                    ylabel(clrbar, 'zscore firing rate');
                    
                    %% plot openfield
%                     sf = subaxis(2,num_trajs+2, [7 8 15 16], 'SpacingVert', Pp.SpVt, ...
%                         'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
%                         Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    % concat, mean, all open field epochs 
%                     xt = openF(ian).outputstack(cl_data_idx(e)).xticks{j};
%                     yt = openF(ian).outputstack(cl_data_idx(e)).yticks{j};
%                     fr{end+1} = [];
%                     fr{end}(xt,yt) = F(ian).outputstack(cl_data_idx(e)).smoothedspikerate{j};
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %d (%s) %d %d %s %s%s',Fp.paths.filenamesave(5:end), animal, ...
                    day, num2str(data_keys(cl_data_idx(:),2)'), ntrode, cluster, Fp.epochEnvironment, filetail, datestr(now, 'yyyymmdd'));
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure(Fp.paths.figdirectory, Fp.paths.filenamesave, sprtit)
                    close all
                end
            end
        end
    end
end


