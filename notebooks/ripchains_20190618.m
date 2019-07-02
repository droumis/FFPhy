

% DO
% swr duration as a function of days
% swr duration as a function of exposure time ( continuous version of the above)
% swr duration as a function of change in learning rate

% load the ripkons for all the animals
% load the behave state for all the animals
% load the linpos for all the animals

%     linpos(ani).animal = animal;
%     linpos(ani).linpos = loaddatastruct(andef{2}, animal, 'linpos');


animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
loadstuff = 1;
plotstuff = 0;

events = struct;
task = struct;
BehaveState = struct;
%% load
if loadstuff
    for ani = 1:length(animals)
        animal = animals{ani};
        andef = animaldef(animal);
        
        % load task, save list of ONLY wtrack days and epochs (not rotated)
        task(ani).animal = animal;
        task(ani).info = loaddatastruct(andef{2}, animal, 'task');
        alleps = evaluatefilter(task(ani).info, 'isequal($environment, ''wtrack'')');
        
        task(ani).wtrackdays = [];
        task(ani).wtrackdayeps = [];
        alldays = unique(alleps(:,1));
        for dy = 1:length(alldays)
            iday_eps = alleps(alleps(:,1)==alldays(dy),:);
            if length(iday_eps(:,1)) == 2 % if there are two wtrack epochs on this day
                task(ani).wtrackdays(end+1) = iday_eps(1,1);
                task(ani).wtrackdayeps = [task(ani).wtrackdayeps; iday_eps];
            end
        end
        
        % load events
        events(ani).animal = animal;
        events(ani).ca1ripples = loaddatastruct(andef{2}, animal, 'ca1rippleskons');
        events(ani).mecripples = loaddatastruct(andef{2}, animal, 'mecrippleskons');
        
        % get rip matrix of [starttime endtime duration day epoch]
        events(ani).ripmatfields = ['start end duration day epoch'];
        events(ani).ca1ripmat = [];
        % ca1 rips
        tmp = arrayfun(@(d,e) [events(ani).ca1ripples{d}{e}{1}.starttime ...
            events(ani).ca1ripples{d}{e}{1}.endtime], task(ani).wtrackdayeps(:,1), ...
            task(ani).wtrackdayeps(:,2), 'un', 0);
        events(ani).ca1ripmat = cell2mat(tmp);
        events(ani).ca1ripmat(:,3) = events(ani).ca1ripmat(:,2)-events(ani).ca1ripmat(:,1);
        events(ani).ca1ripmat(:,[4 5]) = cell2mat(cellfun(@(x,y) repmat(task(ani).wtrackdayeps(y,:), ...
            length(x),1), tmp, num2cell(1:length(tmp))', 'un', 0));
        % mec rips
        events(ani).mecripmat = [];
        tmp = arrayfun(@(d,e) [events(ani).mecripples{d}{e}{1}.starttime ...
            events(ani).mecripples{d}{e}{1}.endtime], task(ani).wtrackdayeps(:,1), ...
            task(ani).wtrackdayeps(:,2), 'un', 0);
        events(ani).mecripmat = cell2mat(tmp);
        events(ani).mecripmat(:,3) = events(ani).mecripmat(:,2)-events(ani).mecripmat(:,1);
        events(ani).mecripmat(:,[4 5]) = cell2mat(cellfun(@(x,y) repmat(task(ani).wtrackdayeps(y,:), ...
            length(x),1), tmp, num2cell(1:length(tmp))', 'un', 0));
        
        % get median, mean, sterr of swr duration per day and per epoch
        for d = 1:length(task(ani).wtrackdays)
            useidx = events(ani).ca1ripmat(:,4) == task(ani).wtrackdays(d);
            durs = events(ani).ca1ripmat(useidx,3);
            events(ani).ca1_durmedian_perday(d) = median(durs);
            events(ani).ca1_durmean_perday(d) = mean(durs);
            events(ani).ca1_durstderr_perday(d) = std(durs) / sqrt(length(durs));

            useidx = events(ani).mecripmat(:,4) == task(ani).wtrackdays(d);
            durs = events(ani).mecripmat(useidx,3);
            events(ani).mec_durmedian_perday(d) = median(durs);
            events(ani).mec_durmean_perday(d) = mean(durs);
            events(ani).mec_durstderr_perday(d) = std(durs) / sqrt(length(durs));
        end

        % load behave state
        tmp = load(sprintf('%s/%sBehaveState.mat',andef{2}, animal));
        BehaveState(ani).animal = animal;
        BehaveState(ani).statechanges = tmp.BehaveState.statechanges;
        BehaveState(ani).statespace = tmp.BehaveState.statespace;
        % get performance (as percent correct) per day and per epoch
        crtde = [BehaveState(ani).statespace.allepsMat(:,7) BehaveState(ani).statespace.allbound(:,[5 6])];
        crtde(~ismember(crtde(:,2), task(ani).wtrackdays),:) = [];
        dayeps = unique(crtde(:,[2 3]), 'rows');
        days = unique(dayeps(:,1));
        
        BehaveState(ani).days = days;
        for d = 1:length(days)
            day = days(d);
            cde_d = crtde(crtde(:,2) == day,:);
            BehaveState(ani).numtrials_perday(d) = length(cde_d(:,1));
            BehaveState(ani).numcorrect_perday(d) = sum(cde_d(:,1));
            BehaveState(ani).perccorrect_perday(d) = sum(cde_d(:,1)) / length(cde_d(:,1));
        end
        
        BehaveState(ani).dayeps = dayeps;
        for de = 1:length(dayeps)
            dayep = dayeps(de,:);
            cde_d = crtde(ismember(crtde(:,[2 3]),dayep,'rows'),:);
            BehaveState(ani).numtrials_perep(de) = length(cde_d(:,1));
            BehaveState(ani).numcorrect_perep(de) = sum(cde_d(:,1));
            BehaveState(ani).perccorrect_perep(de) = sum(cde_d(:,1)) / length(cde_d(:,1));
        end
        % get change in percent correct from each epoch to the next, within
        % each day ( maybe i want to do across day epochs too)
        BehaveState(ani).perccorrectchange_perday = [];
        for d = 1:length(days)
            BehaveState(ani).perccorrectchange_perday(end+1) = diff(...
                BehaveState(ani).perccorrect_perep(dayeps(:,1) == days(d)));
        end
    end
end
%% redo cdf plot of durations over each day for each animal.. 
%%
figure
for ani = 1:length(animals)
    days = unique(events(ani).ca1ripmat(:,4));
    s = subplot(1,length(animals),ani);
    c = jet(length(days));
    for d = 1:length(days)
        ca1dayevents = events(ani).ca1ripmat(events(ani).ca1ripmat(:,4)==days(d),3);
        mecdayevents = events(ani).mecripmat(events(ani).mecripmat(:,4)==days(d),3);        
        ca1dayevents = ca1dayevents(ca1dayevents<.4);
        mecdayevents = mecdayevents(mecdayevents<.4);
%         ca12mecidx = knnsearch(mecdayevents(:,1), ca1dayevents(:,1));
%         ca12mec_dur = ca1dayevents(:,1)-mecdayevents(ca12mecidx,1);
%         v = .4;
%         p = histogram(ca12mec_dur(abs(ca12mec_dur)<v),1000, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
%         p.EdgeColor = c(d,:);
%         cdfplot(ca1dayevents)
%         set(gca,'XScale','log')
%         ecdf(ca1dayevents)
        ecdf(ca1dayevents)
        hold on
        axis tight
        xlabel('swr duration')
%         xlim([-v v])
    end
    title(sprintf('%s',animals{ani}))
    hold off
end
%% detect ripple chains
ca1numchains = nan*ones(length(animals), 10);
mecnumchains = nan*ones(length(animals), 10);
ca1over100 = nan*ones(length(animals), 10);
mecover100 = nan*ones(length(animals), 10);
for ani = 1:length(animals)
   days = unique(events(ani).ca1ripmat(:,4));
   for d = 1:length(days)
        ca1dayevents = events(ani).ca1ripmat(events(ani).ca1ripmat(:,4)==days(d),:);
        ca1dayevents = ca1dayevents(ca1dayevents(:,3)<.4,:);
        mecdayevents = events(ani).mecripmat(events(ani).mecripmat(:,4)==days(d),:);
        mecdayevents = mecdayevents(mecdayevents(:,3)<.4,:);
%         knnsearch(ca1dayevents(:,1), ca1dayevents(:,1));
%         [Idx,D] = rangesearch(ca1dayevents(:,1),ca1dayevents(:,1),.5);
        ca1numchains(ani, d) = length(find(diff(ca1dayevents(:,1)) < .5));
        mecnumchains(ani, d) = length(find(diff(mecdayevents(:,1)) < .5));
        ca1over100(ani,d) = length(find(ca1dayevents(:,3)>.100)) / length(ca1dayevents(:,1));
        mecover100(ani,d) = length(find(mecdayevents(:,3)>.100)) / length(mecdayevents(:,1));
%         for i = 1:length(Idx)
%             if ismember(Idx{i}, Idx{i+1:end})
%         end
%         length(cellfun(@length, Idx) > 1)
    end
end
figure
plot(mecover100')
%% Learning rate
for ani = 1:length(animals)
   BehaveState(ani).statespace
    
end





%%
%
% %% get the rip duration distribution for each day
% % day epoch startrip endrip
% for ani = 1:length(animals)
%     wtrackeps = evaluatefilter(task(ani).info, 'isequal($environment, ''wtrack'')');
%     days = unique(wtrackeps(:,1));
%     task(ani).wtrackdays = [];
%     for dy = 1:length(days)
%         numweps = length(find(wtrackeps(:,1)==days(dy)));
%         if numweps == 2
%             task(ani).wtrackdays(end+1) = days(dy);
%         end
%     end
%     task(ani).wtrackeps = wtrackeps;
%     events(ani).wtrackeps = wtrackeps;
%     events(ani).wtrackrips = arrayfun(@(d,e) [events(ani).ca1ripples{d}{e}{1}.starttime events(ani).ca1ripples{d}{e}{1}.endtime], wtrackeps(:,1), wtrackeps(:,2), 'un', 0);
%     events(ani).wtrackripdur = arrayfun(@(d,e) [events(ani).ca1ripples{d}{e}{1}.endtime - events(ani).ca1ripples{d}{e}{1}.starttime], wtrackeps(:,1), wtrackeps(:,2), 'un', 0);
%     events(ani).wtrackripdur_perday = {};
%     events(ani).wtrackrips_perday = {};
%     for d = 1:length(task(ani).wtrackdays)
%         day = task(ani).wtrackdays(d);
%         eps2dayidx = find(events(ani).wtrackeps(:,1) == day);
%         events(ani).wtrackripdur_perday{d} = cell2mat(events(ani).wtrackripdur(eps2dayidx));
%         events(ani).wtrackrips_perday{d} = cell2mat(events(ani).wtrackrips(eps2dayidx));
%     end
% end
% %% get mean, sterr of duration per day
% for ani = 1:length(animals)
%     usdesidx = find(ismember(events(ani).wtrackeps(:,1), task(ani).wtrackdays));
%     events(ani).durmedian_perep = cell2mat(cellfun(@median,events(ani).wtrackripdur(usdesidx), 'un', 0));
%     events(ani).durmean_perep = cell2mat(cellfun(@mean,events(ani).wtrackripdur(usdesidx), 'un', 0));
%     events(ani).durstd_perep = cell2mat(cellfun(@std,events(ani).wtrackripdur(usdesidx), 'un', 0));
%     events(ani).stderror_perep = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)), events(ani).wtrackripdur(usdesidx), 'un', 0));
%
%     events(ani).durmedian_perday = cell2mat(cellfun(@median,events(ani).wtrackripdur_perday(task(ani).wtrackdays), 'un', 0));
%     events(ani).durmean_perday = cell2mat(cellfun(@mean,events(ani).wtrackripdur(task(ani).wtrackdays), 'un', 0));
%     events(ani).durstd_perday = cell2mat(cellfun(@std,events(ani).wtrackripdur(task(ani).wtrackdays), 'un', 0));
%     events(ani).stderror_perday = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)), events(ani).wtrackripdur(task(ani).wtrackdays), 'un', 0));
% end
if plotstuff
    %% plot median duration per ep (change this to per day)
    ifig = figure('units','normalized','position',[.1 .1 .8 .5]);
    for ani = 1:length(animals)
        subplot(2,ceil(length(animals)/2),ani);
        e = errorbar(1:length(events(ani).ca1_durmedian_perday), events(ani).ca1_durmedian_perday,events(ani).ca1_durstderr_perday);
        e.LineWidth = 1.5;
        e.Color = 'k';
        ylabel('duration s')
        xlabel('day')
        axis tight
        xticks(1:length(events(ani).ca1_durmedian_perday))
        title(animals(ani))
    end
    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('median ca1 swr duration');
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center','FontSize', 12);
    
        %% plot median duration per day
    ifig = figure('units','normalized','position',[.1 .1 .5 .5]);
    r = nan*ones(length(animals),8);
    for ani = 1:length(animals)
%       subplot(1,2,1);
      r(ani,1:length(events(ani).ca1_durmedian_perday)) = events(ani).ca1_durmedian_perday;
        e = errorbar(1:length(events(ani).ca1_durmedian_perday), events(ani).ca1_durmedian_perday,events(ani).ca1_durstderr_perday);
        hold on
        e.LineWidth = 1.5;
%         e.Color = 'k';
        ylabel('SWR duration s')
        xlabel('day')
        axis tight
%         xticks(1:length(events(ani).ca1_durmean_perday))
%         title(animals(ani))
    end
    legend(animals{:})
    ylim([0.02 .15])
%     
%     subplot(1,2,2);
%     mall = nanmean(r,1)
%     plot(mall(1:5))

    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('median ca1 swr duration');
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center','FontSize', 12);
    
            %% plot median MEC SWR duration per day
    ifig = figure('units','normalized','position',[.1 .1 .5 .5]);
    r = nan*ones(length(animals),8);
    for ani = 1:length(animals)
%       subplot(1,2,1);
      r(ani,1:length(events(ani).mec_durmedian_perday)) = events(ani).mec_durmedian_perday;
        e = errorbar(1:length(events(ani).mec_durmedian_perday), events(ani).mec_durmedian_perday,events(ani).mec_durstderr_perday);
        hold on
        e.LineWidth = 1.5;
%         e.Color = 'k';
        ylabel('SWR duration s')
        xlabel('day')
        axis tight
%         xticks(1:length(events(ani).ca1_durmean_perday))
%         title(animals(ani))
    end
    legend(animals{:})
    ylim([0.02 .15])
%     
%     subplot(1,2,2);
%     mall = nanmean(r,1)
%     plot(mall(1:5))

    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('median mec swr duration');
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center','FontSize', 12);
    %% Pool the median SWR durations across animals, then test significance using kruskalwallis
    ca1use = {};
    mecuse = {};
    for ani = 1:length(animals)
        ca1use{ani} = events(ani).ca1ripmat(events(ani).ca1ripmat(:,3) < 400,:);
        mecuse{ani} = events(ani).mecripmat(events(ani).mecripmat(:,3) < 400,:);
    end
    
    ca1usemat = cell2mat(ca1use');
    mecusemat = cell2mat(mecuse');
    %%
    andymedca1 = nan*ones(length(animals), 5);
    andysterca1 = nan*ones(length(animals), 5);
    andymedmec = nan*ones(length(animals), 5);
    andystermec = nan*ones(length(animals), 5);
    for d = 1:length(1:5)
        c = ca1usemat(ca1usemat(:,4)==d,3);
        m = mecusemat(mecusemat(:,4)==d,3);
        ca1medpd(d) = median(c);
        ca1sterpd(d) = std(c) / sqrt(length(c));
        mecmedpd(d) = median(m);
        mecsterpd(d) = std(m) / sqrt(length(m));
        
        for ani = 1:length(animals)
           adca1 = ca1use{ani}(ca1use{ani}(:,4)==d,3);
           admec = mecuse{ani}(mecuse{ani}(:,4)==d,3);
           andymedca1(ani,d) = median(adca1);
           andysterca1(ani, d) = std(adca1) / sqrt(length(adca1));
           andymedmec(ani,d) = median(admec);
           andystermec(ani, d) = std(admec) / sqrt(length(admec));
        end
    end
    %%
        figure
    errorbar(mean(andymedca1), std(andymedca1) / sqrt(length(andymedca1)), 'LineWidth', 1.2, 'Color', 'k')
    xlim([.5 5.5])
    title('CA1 SWR duration mean of animal median per day, stderr')
    xlabel('day')
    ylabel('duration s')
    ylim([0.02 .15])
    
    figure
    errorbar(mean(andymedmec), std(andymedmec) / sqrt(length(andymedmec)), 'LineWidth', 1.2, 'Color', 'k')
    xlim([.5 5.5])
    title('MEC SWR duration mean of animal median per day, stderr')
    xlabel('day')
    ylabel('duration s')
    ylim([0.02 .15])
    
    
    
    figure
    errorbar(ca1medpd, ca1sterpd, 'LineWidth', 1.2, 'Color', 'k')
    xlim([.5 5.5])
    title('CA1 SWR duration median, stderr')
    xlabel('day')
    ylabel('duration s')
    
    figure
    errorbar(mecmedpd, mecsterpd, 'LineWidth', 1.2, 'Color', 'k')
    xlim([.5 5.5])
    title('MEC SWR duration median, stderr')
    xlabel('day')
    ylabel('duration s')
    
%     pc = kruskalwallis(ca1medpd, );
%     pm = kruskalwallis(mecmedpd, mecsterpd);
    pc = kruskalwallis(ca1usemat(:,3), ca1usemat(:,4));
    pm = kruskalwallis(mecusemat(:,3), mecusemat(:,4));
%     
    pc = kruskalwallis(andymedca1);
    pm = kruskalwallis(andymedmec);
    
    
    %% plot the entire distribution for each day per animal
    
    %% plot performance per day
    
    for ani = 1:length(animals)
        subplot(2,ceil(length(animals)/2),ani)
        usdys = find(ismember(BehaveState(ani).days, task(ani).wtrackdays));
        days = BehaveState(ani).days(usdys);
        plot(days, BehaveState(ani).perccorrect_perday(usdys))
        hold on
        title(animals{ani})
        xlabel('wtrackday')
        ylabel('% correct')
    end
    %% plot performance per ep
    for ani = 1:length(animals)
        subplot(2,ceil(length(animals)/2),ani)
        usdys = find(ismember(BehaveState(ani).dayeps(:,1), task(ani).wtrackdays));
        dayeps = BehaveState(ani).dayeps(usdys,:);
        plot((.5:.5:length(dayeps)/2), BehaveState(ani).perccorrect_perep(usdys))
        title(animals{ani})
        xlabel('wtrackdayep')
        ylabel('% correct')
        hold off
    end
    %% plot the change in perc correct between epochs within each day
    ifig2 = figure('units','normalized','position',[.1 .1 .8 .5]);
    for ani = 1:length(animals)
        subplot(2,ceil(length(animals)/2),ani)
        usdyeps = find(ismember(BehaveState(ani).dayeps(:,1), task(ani).wtrackdays));
        dayeps = BehaveState(ani).dayeps(usdyeps,:);
        dyas = unique(dayeps(:,1));
        BehaveState(ani).perccorrectchange_perday = [];
        for dy = 1:length(dyas)
            BehaveState(ani).perccorrectchange_perday(end+1) = diff(...
                BehaveState(ani).perccorrect_perep(dayeps(:,1) == dyas(dy)));
        end
        p = plot(dyas, BehaveState(ani).perccorrectchange_perday);
        p.Color = 'b';
        p.LineWidth = 1.5;
        title(animals{ani})
        xlabel('wtrackday')
        ylabel('change % correct')
        axis tight
        xticks(dyas)
    end
    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig2);
    sprtit = sprintf('percent correct change within day');
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center','FontSize', 12);
    %% plot regression between median swr duration and perf change per day
    ifig3 = figure('units','normalized','position',[.1 .1 .8 .5]);
    for ani = 1:length(animals)
        subplot(2,ceil(length(animals)/2),ani)
        Y = BehaveState(ani).perccorrectchange_perday;
        X = events(ani).ca1_durmedian_perday;
        scatter(X,Y,30, copper(length(Y)),'filled')
        p = polyfit(X,Y,1);
        f = polyval(p,X);
        hold on
        plot(X,f,'-b')
        
        title(animals{ani})
        xlabel('% correct change')
        ylabel('SWR duration')
        axis tight
    end
    
    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig3);
    sprtit = sprintf('ca1 perccorrect change X swr duration');
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center','FontSize', 12);
end
%% 

%% Get iterval between the start of each CA1 SWR to the nearest MEC SWR
% then group into per epoch or per day distribution and plot
% then regression against percent correct, percent correct change
% time since start of day? 1st epoch vs 2nd epoch essentially
% correct vs incorrect trial ( this is a bit weird bc post-rewarded swrs
% would be assigned to the category that follows rather than the trial leading up to it)
% time since last well visit? or maybe just at well vs away from well

for ani = 1:length(animals)
    days = unique(events(ani).ca1ripmat(:,4));
    s = subplot(1,length(animals),ani);
    c = jet(length(days));
    for d = 1:length(days)
        ca1dayevents = events(ani).ca1ripmat(events(ani).ca1ripmat(:,4)==days(d),:);
        mecdayevents = events(ani).mecripmat(events(ani).mecripmat(:,4)==days(d),:);
        ca12mecidx = knnsearch(mecdayevents(:,1), ca1dayevents(:,1));
        ca12mec_dur = ca1dayevents(:,1)-mecdayevents(ca12mecidx,1);
        v = .4;
        p = histogram(ca12mec_dur(abs(ca12mec_dur)<v),1000, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
        p.EdgeColor = c(d,:);
        hold on
%         p.EdgeAlpha = 0;
%         p.EdgeAlpha = 0;
        xlim([-v v])
%         ylim([0 20])
    end
    title(sprintf('%s',animals{ani}))
    hold off
end
%% test whether the duration of individual ripples correlates with interval
% to closest mec ripple
for ani = 1:length(animals)
    days = unique(events(ani).ca1ripmat(:,4));
    s = subplot(1,length(animals),ani);
    c = jet(length(days));
    for d = 1:length(days)
        ca1dayevents = events(ani).ca1ripmat(events(ani).ca1ripmat(:,4)==days(d),:);
        ca1dayevents = ca1dayevents(ca1dayevents(:,3) < .25,:);
        mecdayevents = events(ani).mecripmat(events(ani).mecripmat(:,4)==days(d),:);
        mecdayevents = mecdayevents(mecdayevents(:,3) < .25,:);
        ca12mecidx = knnsearch(mecdayevents(:,1), ca1dayevents(:,1));
        ca12mec_dur = abs(ca1dayevents(:,1)-mecdayevents(ca12mecidx,1));
        ca12mec_dur = ca12mec_dur(ca12mec_dur < 10,:);
        ca1dayevents = ca1dayevents(ca12mec_dur < 10,:);
        scatter(ca12mec_dur, ca1dayevents(:,3), '.')
        xlabel('dur ca1 2 mec')
        ylabel('ca1 dur')
%         v = .5;
%         p = histogram(ca12mec_dur(abs(ca12mec_dur)<v),1000, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
%         p.EdgeColor = c(d,:);
        hold on
%         p.EdgeAlpha = 0;
%         p.EdgeAlpha = 0;
%         xlim([-v v])
%         ylim([0 20])
    end
    title(sprintf('%s',animals{ani}))
    hold off
end

%% test whether ripple chains (how to detect these??) translate into the likelihood of being near mec rip (i.e. propagation)

for ani = 1:length(animals)
    days = unique(events(ani).ca1ripmat(:,4));
    iso = [];
    cha = [];
    subplot(2,ceil(length(animals)/2),ani)
    for d = 1:length(days)
        ca1dayevents = events(ani).ca1ripmat(events(ani).ca1ripmat(:,4)==days(d),:);
        a = [inf; diff(ca1dayevents(:,1))];
        chains = ca1dayevents(a < .25,:);
        chains2mecidx = knnsearch(mecdayevents(:,1), chains(:,1));`
        chains2mec_dur = chains(:,1)-mecdayevents(chains2mecidx,1);
        
        isos = ca1dayevents(a(2:end) > .25,:);
        iso2mecidx = knnsearch(mecdayevents(:,1), isos(:,1));
        iso2mec_dur = isos(:,1)-mecdayevents(iso2mecidx,1);

        iso = [iso;iso2mec_dur];
        cha = [cha;chains2mec_dur];``
%         figure
%         histogram(chains2mec_dur,100)
%         [h,p,ci,stats] = ttest2(chains2mec_dur, iso2mec_dur);
%         ca1dayevents = ca1dayevents(ca1dayevents(:,3) < .25,:);
    end
    boxplot([iso;cha], [zeros(length(iso),1); ones(length(cha),1)])
end
%% cumulative reward history
figure
for ani = 1:length(animals)
    set(gca, 'color', 'white')
%     f1 = subplot(1,2,1);
    semilogy(cumsum(BehaveState(ani).statespace.allepsMat(:,7)), 'LineWidth', 1.2)
    xlim([0 100])
    ylim([0 100])
    xlabel('well visit #')
    ylabel('log cumulative reward')
    hold on
%     f2 = subplot(1,2,2);
%     loglog(cumsum(BehaveState(ani).statespace.allepsMat(:,7)), 'LineWidth', 1.2)
%     xlabel('log trial #')
%     ylabel('log cumulative reward')
%     hold on
end
hold off
legend(animals{:})












