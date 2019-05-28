

% get ripple locked ripple triggered LFP phase for each ntrode

animals = {'D10'};

filtfunction = 'riptriglfp';
env = 'wtrack';

runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
combineData = 0; 
% only using JZ1 < day 8 bc day 8 missing data issue still open
% JZ3 day 3 is also missing

plotfigs = 0;
pausefigs = 0;
savefigs = 1;


%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework == 1
    paths = make_paths(Fp.filtfunction, Fp.epochEnvironment);
    Fp = load_filter_params({env, 'ripples', filtfunction, '2s_win'});
    Fp.animals = animals;
    F = createfilter('animal', Fp.animals, Fp.createfilterset{:});
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    for a = 1:length(F) % save filter detailes along with results
        F(a).datafilter_params = Fp;
    end
    F = runfilter(F);
end
%% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput == 1
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave)
end
%% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput == 1
    Fp = load_filter_params({env, 'ripples', filtfunction});
    Fp.animals = animals;
    paths = make_paths(Fp.filtfunction, Fp.epochEnvironment);
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, Fp.animals);
end
%%
% ok so it looks the results are saved into a output/per-day cell array,
% this contains a struct array with epochs for rows in there ep ID row #
% it contains fields: data, which is a per-eeg type cell array,
% that contains a per-ripple cell array

%I think maybe instead I just want to rely on the day/epoch keys i would generate from
%either the task/cell/tetinfo or from the data itself..

eegtype = 1;
if combineData
    lfp = struct;
    for ian = 1:length(F)
        lfp(ian).animal = F{ian}.F.animal{3};
        lfp(ian).lfptypes = F{ian}.F.datafilter_params.LFPtypes;
        lfp(ian).LFPrangesHz = F{ian}.F.datafilter_params.LFPrangesHz;
        lfp(ian).data_dims = {'swr#', 'time', 'ntrode'};
        lfp(ian).dayeps = F{ian}.F.epochs{1, 1};
        if lfp(ian).animal == 'JZ4'
            lfp(ian).dayeps = lfp(ian).dayeps(lfp(ian).dayeps(:,1)<8, :);
        elseif lfp(ian).animal == 'JZ3'
            lfp(ian).dayeps = lfp(ian).dayeps(lfp(ian).dayeps(:,1)<3, :);
        end
        lfp(ian).data = {}; lfp(ian).numrips = [];
        for de = 1:length(lfp(ian).dayeps(:,1)) % day epoch
            day = lfp(ian).dayeps(de,1);
            epoch = lfp(ian).dayeps(de,2);
            for t = 1:length(lfp(ian).lfptypes) % LFP type
                try
                % rip# X time X ntrode
                tmp = F{ian}.F.output{day}(epoch).data{t};
                lfp(ian).data{t}{de} = permute(cat(3,tmp{:}), [3 2 1]);
                catch
                    continue
                end
            end
            lfp(ian).numrips{de} = length(lfp(ian).data{t}{de}(:,1,1));
        end
        for t = 1:length(lfp(ian).lfptypes) % LFP type
            lfp(ian).data{t} = cell2mat(lfp(ian).data{t}'); % stack across time
        end
        ntrodes = cellfun(@(x) x(:,3), {F{ian}.F.output{day}(epoch).index}, 'un', 0);
        lfp(ian).ntrodes = ntrodes{1};
        lfp(ian).numrips = cell2mat(lfp(ian).numrips);
    end
end
%% PLOT
if plotfigs
    Pp = load_plotting_params({'riptriglfp_acrossdays'});
    for ian = 1:length(lfp)
        for nt = 1:length(lfp(ian).ntrodes)
            if savefigs && ~pausefigs
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                clf
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            animal = lfp(ian).animal;
            ntrode = lfp(ian).ntrodes(nt);
            try
                for t = 1:length(lfp(ian).lfptypes)
                    
                    subaxis(1,length(lfp(ian).lfptypes), t);
                    %                 indwin = plotwin*1500;
                    a = lfp(ian).data{t}(:,:,ntrode);
                    mid = ceil(size(a,2)/2);
                    pwin = [.5 .5];
                    srate = 1500;
                    a = double(a(:,mid-(pwin(1)*srate):mid+(pwin(2)*srate)));
                    m = nanmean(a,2);
                    s = nanstd(a, [], 2);
                    imagesc((a-m)./s)
                    %                 imagesc(lfp(ian).data{t}(:,:,ntrode));
                    % day/epoch lines
                    dayind = find(diff(lfp(ian).dayeps(:,1))>0);
                    epbounds = cumsum(lfp(ian).numrips);
                    daybounds = epbounds(dayind);
                    line([1 size(lfp(ian).data{t},2)], [epbounds; epbounds], 'color', [.8 .8 .8])
                    line([1 size(lfp(ian).data{t},2)], [daybounds; daybounds], 'color', [0 0 0])
                    line([ceil(mid/2) ceil(mid/2)], [1 size(lfp(ian).data{t},1)], 'color', [0 0 0], 'linestyle', '--')
                    title(lfp(ian).lfptypes{t})
                    if t == 1
                        xticks([1 ceil(mid/2) mid])
                        xticklabels({'-.5','0','.5'})
                        ylabel('swr#')
                    end
                end
                cmap = colormap('parula');
                colormap(cmap(end:-1:1,:));
                %% ---- super title and colorbar----
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('riptriglfp zscore alldays %s %s nt%d', Fp.epochEnvironment, ...
                    animal, ntrode);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', ...
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
            catch
                continue
                fprintf('error trying to plot %s nt%d.. skipping \n', animal, ntrode);
            end
        end
    end
end








