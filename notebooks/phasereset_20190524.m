
%{

Jutras Fries Bufallo 2013
- Powervalues were calculated for the presaccade period (the 600-msperiod preceding
saccade onset) and the postsaccade period (the600-ms period following a
400-ms postsaccade“buffer”period)using separate methods to obtain individual trial power
measuresand trial-averaged power measures.

For individual trial power measures:
 - multplied the raw LFP traces for each 600 ms period with one orthogonal
 taper function before Fourier transofrmation, providing spectral smoothing
 of +-1.67Hz.
- We then calculated the average power in the 6.7 to 11.6 Hz freq band,
across trials, to obtain a theta power value for each LFP for presaccade
and postsaccade periods.

For trial averaged power measures, the same method was yused, with the
exception that before taper multipplication and Fourier transformation,
we calculated the average LFP signal (i.e. evoked signal) across trials.
-.. then applied taper multiplication and fouruer transformation to each
LFP segment (trial averaged pre and post saccade) and caluclated the
average power in the 6.7 to 11.6 Hz freq band, producing a theta power
value for pre and post trial averaged signals.

A paired t test was used to determine whether theta powrer was
significantly different for pre and post periods, serperately for
individual trial theta power values and trial avergaed theta power values.
%}

% 1. load the ripple triggered LFP
% 2. collect results per ntrode
% - what are my other functions i've made to do something similar, either
% - combine by days, epochs, ntrodes, epochtype_perday, etc..
% - i need to centralize these whether it's with spikes or lfp or
% whatever.. which means standardizing the output of the analysis
% functions..
% 3. average LFP across rips
% 4. calculate theta power pre and post swr

% get ripple locked ripple triggered LFP phase for each ntrode

animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};

filtfunction = 'riptriglfp';
env = 'wtrack';

loadFilterOutput = 0;
combineData_perNtrode = 0;
% only using JZ4 < day 8 bc day 8 missing data issue still open
% JZ3 day 3 is also missing

plotfigs = 0;
pausefigs = 0;
savefigs = 1;

%% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput == 1
    Fp = load_filter_params({env, 'ripples', filtfunction});
    Fp.animals = animals;
    paths = make_paths(Fp.filtfunction, Fp.epochEnvironment);
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, Fp.animals);
end
%% ---- compile results from {day}(ep).data(<ntXv>) into lfp(animal).data(<ripXtimeXnt>)
% with day and ep info along at the level under ntrode.. i.e.
% ntrode.dayeps()
eegtype = 1;
if combineData_perNtrode
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
        lfp(ian).data = {}; lfp(ian).numrips = []; lfp(ian).time = [];
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
            % collect day/epoch level info
            lfp(ian).numrips{de} = length(lfp(ian).data{t}{de}(:,1,1));
            lfp(ian).ripDayEp{de} = repmat([day epoch], ...
                length(lfp(ian).data{t}{de}(:,1,1)), 1);
            lfp(ian).ripStartIdx{de} = F{ian}.F.output{day}(epoch).eventStartIndices;
            lfp(ian).ripEndIdx{de} = F{ian}.F.output{day}(epoch).eventEndIndices;
            
            lfp(ian).ripStartTime{de} = F{ian}.F.output{day}(epoch).LFPtimes(...
                F{ian}.F.output{day}(epoch).eventStartIndices);
            lfp(ian).ripEndTime{de} = F{ian}.F.output{day}(epoch).LFPtimes(...
                F{ian}.F.output{day}(epoch).eventEndIndices);
            
        end
        for t = 1:length(lfp(ian).lfptypes) % LFP type
            lfp(ian).data{t} = cell2mat(lfp(ian).data{t}'); % stack across time
        end
        ntrodes = cellfun(@(x) x(:,3), {F{ian}.F.output{day}(epoch).index}, 'un', 0);
        lfp(ian).ntrodes = ntrodes{1};
        lfp(ian).numrips = cell2mat(lfp(ian).numrips);
    end
end


%% PLIT
if plotfigs
    
    %% create a morlet wavelet
    
    srate = 1500;
    time = -1:1/srate:1;
    frex = 7;
    % complex sin wave
    sine_wave = exp( 1i*2*pi*frex.*time );
    % gaussian window
    numcycles = 7;
    win = numcycles / (2*pi*frex);
    gaus_win = exp( (-time.^2) ./ (2*win*2) );
    % complex morlet wavelet
    cmw = sine_wave .* gaus_win;
    
    %%
    % riptriglfp data are in ripnum, time, ntrode dims
    nt = 11;
    exampleripnum = 102;
    ian = 2;
    
    % gather data for this ntrode/ripnum
    exampletrace = lfp(ian).data{1}(exampleripnum,:,nt);
    ripDayEps = cell2mat(lfp(ian).ripDayEp');
    iripDayEp = ripDayEps(exampleripnum,:);
    ripStarts = cell2mat(lfp(ian).ripStartTime');
    ripEnds = cell2mat(lfp(ian).ripEndTime');
    iripStart = ripStarts(exampleripnum);
    iripEnd = ripEnds(exampleripnum);
    
    win = F{ian}.F.output{day}(epoch).win; % hack cuz the win does not change
    time = -win(1):1/srate:win(2);
    
    plot(time, exampletrace)
    xlabel('time (s)')
    ylabel('uV')
    title(sprintf('%s %d %d nt%d rip%d %.03f-%.03f s', animals{ian}, iripDayEp(1), ...
        iripDayEp(2), nt, exampleripnum, iripStart, iripEnd))
    %% plot rip mean traces
    % i want to plot patches of all the ripples within this window
    % look at how i did it with riptriglfp bc i want to also include the
    % peak std in the plot
    win = F{ian}.F.output{1}(2).win; % hack cuz the win does not change
    time = -win(1):1/srate:win(2);
    for ian = 1:length(animals)
        figure(ian)
        % gather data for this ntrode/ripnum
        meanXrips = squeeze(mean(lfp(ian).data{1}(:,:,:), 1))';
        waterfall(time, lfp(ian).ntrodes, meanXrips)
        %     plot(time, )
        xlabel('time s')
        ylabel('ntrode')
        zlabel('LFP uV')
        sprtit = sprintf('%s mean riptriglfp alleps 3d wtrack', animals{ian});
        title(sprtit)
        save_figure(paths.figdirectory, paths.filenamesave, sprtit)
        pause
    end
    
    % so it looks like on average.. there is theta power both before and
    % after.. i need to look at this over ripples over time.. see what the
    % right thing to do is to characterize the shift in the phase reset
    
    % also there's def something wrong with D12...
    
    %%
    %     for ian = 1:length(animals)
    lfptypes = F{ian}.F.output{1}(2).LFPtypes;
    for f = 1:5
    for ian = 1:length(animals)
        Pp = load_plotting_params({'riptriglfp_acrossdays'});
        if savefigs && ~pausefigs
            close all
            ifig = figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            clf
            ifig = figure('units','normalized','position',Pp.position);
        end
        set(gcf,'color','white')
        for n = 1:15
%             figure(1)
%             subaxis(15,1,n,'Spacing', .01, 'Padding', 0, 'Margin', .03)
%             ntrips = lfp(ian).data{1}(:,:,n);
%             image(ntrips)
%             ylabel(sprintf('%d',n))
%             set(gca, 'xTick', [])
%             set(gca, 'YTick', [])
%             figure(2)
%             subaxis(1,15,n,'Spacing', .01, 'Padding', 0, 'Margin', .03)
%             ntrips = lfp(ian).data{1}(:,:,n);
%             image(ntrips)
%             xlabel(sprintf('%d',n))
%             set(gca, 'xTick', [])
%             set(gca, 'YTick', [])
%             figure(3)
            subaxis(4,4,n,'Spacing', .02, 'Padding', .01, 'Margin', .04)
            ntrips = lfp(ian).data{f}(:,:,n);
            image(ntrips)
            ylabel(sprintf('%d',n))
            set(gca, 'xTick', [])
            set(gca, 'YTick', [])
        end
        cmap = colormap('bone');
        colormap(cmap(end:-1:1,:));
        % ---- super title and colorbar----
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        sprtit = sprintf('%s riptriglfp %sband %s alleps allmec',...
            animals{ian}, lfptypes{f}, Fp.epochEnvironment);
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center');
        
        % ---- pause, save figs ----
        if pausefigs
            pause
        end
        if savefigs
            save_figure(paths.figdirectory, paths.filenamesave, sprtit)
        end
    end
    end
    
    %%
end


















