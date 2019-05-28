

%{
- calculate mu corr measures per rip

%}

Fp = load_filter_params({'wtrack','ripples','dfa_perripspikingcorr'});
Fp.animals = {'D12'}; %D12
Fp.days = [1:6]; % don't use D12 d7 until sorting is done
% Fp.days = [2];
env = 'wtrack';

runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;

plotstuff = 0;
savefigs = 0;
pausefigs = 0;

%% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp.filtfunction, Fp.epochEnvironment);
%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework == 1
    F = createfilter('animal', Fp.animals, 'epochs', ...
        Fp.epochfilter, 'tetrodepairs', Fp.tetpairfilter, 'iterator', ...
        Fp.iterator, 'excludetimefilter', Fp.timefilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes);
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
    for ian = 1:length(Fp.animals)
        F(ian) = load_filter_output(paths.filtOutputDirectory, ...
            paths.filenamesave, 'filetail', sprintf('_%s.mat',Fp.animals{ian}));
    end
end

if plotstuff
    
    Pp = load_plotting_params('ripcorr_Xdays');
    % so now that i have a few measures of rip trig spike correlation, i
    % want to collect each into a distribution per day per pair.
    % probably the best way to do this is a box whisker or violin
    animals = cellfun(@(x) x{1}, {F.animal}', 'un', 0);
    for an = 1:length(animals)
        data_indices = cell2mat({F(an).output{1}.index}');
        % get mapping to each unique tetrode
        [pairs,~, data2pairs] = unique(data_indices(:,[3 4]), 'rows');
        for ipair = 1:length(pairs(:,1))
            %init fig
            if savefigs && ~pausefigs;
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                clf
                ifig = figure(1);%'units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            % collect data for this pair
            
            pair_datainds = find(data2pairs == ipair);
            ntpair_ids = data_indices(pair_datainds,:);
            days = unique(ntpair_ids(:,1));
            ec = []; raw_corr = []; instaFR=[]; bin_10ms=[]; rms=[];
            for	d = 1:length(days)
                % combine across epochs for this day
                % find this day for tet inds into data
                ipairday_ids = ntpair_ids(find(ntpair_ids(:,1)==days(d)),:);
                ipair_data_inds = find(ismember(data_indices, ...
                    ipairday_ids,'rows'));
                for ipair_day_ep_ind = 1:length(ipair_data_inds(:,1))
                    ide = F(an).output{1}(ipair_data_inds(ipair_day_ep_ind));
                    daygrp = days(d)*ones(length(ide.ntAB_excesscorr),1);
                    try
                        ec = [ec; [ide.ntAB_excesscorr' daygrp]];
                    catch
                        ec = [ec; [ide.ntAB_excesscorr daygrp]];
                    end
                    try
                        raw_corr = [raw_corr; [ide.ntAB_spike_corr daygrp]];
                        instaFR = [instaFR; [ide.ntAB_instaFR_corr daygrp]];
                        bin_10ms = [bin_10ms; [ide.ntAB_10msFR_corr daygrp]];
                        rms = [rms; [ide.ntAB_rms' daygrp]];
                    catch
                        raw_corr = [raw_corr; [ide.ntAB_spike_corr' daygrp]];
                        instaFR = [instaFR; [ide.ntAB_instaFR_corr' daygrp]];
                        bin_10ms = [bin_10ms; [ide.ntAB_10msFR_corr' daygrp]];
                        rms = [rms; [ide.ntAB_rms' daygrp]];
                    end
                end
%                 ech = subaxis(5,1,1); 
% %                 a = histogram(ec(ec(:,2)==days(d),1), 100);
% %                 a(1).FaceAlpha = .5;
% %                 subaxis(5,1,2)
%                 p = plot(ec(ec(:,2)==days(d),1));
%                 p.Color = Pp.cmap(d,:);
%                 hold on
            end
%             hold off
%             axis tight
            ech = subaxis(5,1,1); 
            boxplot(ec(:,1), ec(:,2), 'Notch','on')
            ylabel('excess corr')
            rch = subaxis(5,1,2);
            boxplot(raw_corr(:,1), raw_corr(:,2)) % second arg is day grouping
            ylabel('corr')
            insh = subaxis(5,1,3);
            boxplot(instaFR(:,1), instaFR(:,2)) % second arg is day grouping
            ylabel('instaFR corr')
            b10h = subaxis(5,1,4);
            boxplot(bin_10ms(:,1),bin_10ms(:,2)) % second arg is day grouping
            ylabel('10msbin corr')
            rmh = subaxis(5,1,5);
            boxplot(rms(:,1), rms(:,2)) % second arg is day grouping
            ylabel('rms')
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('riptrig MU corr %s %s - nt%d nt%d', env, ...
                animals{an}, pairs(ipair,1), pairs(ipair,2));
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
        end
    end
end




















