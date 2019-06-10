

animals = {'JZ1'};
add_params = {'wtrack'};

run_lfp_ff = 0;
save_lfp = run_lfp_ff;

run_spikes_ff = 0;
save_spikes = run_spikes_ff;

loadstuff = 1;
load_lfp = loadstuff;
load_spikes = loadstuff;
load_pos = 1;
load_rips = 1;

stackstuff = 1;
stack_lfp = stackstuff;
stack_spikes = stackstuff;

plotfigs = 0;
pausefigs = 0;
savefigs = 0;

%% run lfp filter/func
Fp.animals = animals;
Fp.add_params = add_params;
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);

if run_lfp_ff == 1
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F);
    for d = 1:length(F)
        F(d).datafilter_params = Fp;
    end
end
%% save lfp data
if save_lfp == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        ['_', Fp.epochEnvironment])
end
%% load lfp
if load_lfp
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_', Fp.epochEnvironment]);
end
%% stack LFP
if stack_lfp
    lfpstack = stack_riptriglfp(F);
end

%% run spiking FF
sFp.animals = animals;
sFp.add_params = add_params;
sFp.filtfunction = 'dfa_riptrigspiking';
sFp = load_filter_params(sFp, 'add_params',sFp.add_params);

if run_spikes_ff == 1
    spikesF = createfilter('animal',sFp.animals, 'epochs', sFp.epochfilter, ...
        'cells',sFp.cellfilter, 'excludetime',sFp.timefilter, 'iterator',sFp.iterator);
    spikesF = setfilterfunction(spikesF,sFp.filtfunction,sFp.datatypes,sFp.options{:});
    spikesF = runfilter(spikesF);
    for d = 1:length(spikesF)
        spikesF(d).datafilter_params = sFp;
    end
end
%% save spikes data
if save_spikes == 1
    save_data(spikesF,sFp.paths.filtOutputDirectory,sFp.paths.filenamesave, 'filetail',...
        ['_',sFp.epochEnvironment])
end
%% load spikes
if load_spikes
    spikesF = load_filter_output(sFp.paths.filtOutputDirectory,sFp.paths.filenamesave, ...
        sFp.animals, 'filetail', ['_',sFp.epochEnvironment]);
end

%% stack spikes
if stack_spikes
    spikestack = stack_riptrigspiking(spikesF);
end

%% load pos
if load_pos
    for ian = 1:length(Fp.animals)
        andef = animaldef(Fp.animals{ian});
        pos{ian} = loaddatastruct(andef{2},andef{3},'pos');
        linpos{ian} = loaddatastruct(andef{2},andef{3},'linpos');
        behave{ian} = load(sprintf('%s/%sBehaveState.mat',andef{2}, andef{3}));
    end
end

%% load events
if load_rips
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        %         ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        %         ntinfoAll = cellfetch(ntinfo, '', 'alltags', 1);
        rips{ian} = loaddatastruct(aninfo{2}, animal, 'ca1rippleskons');
        noise{ian} = loaddatastruct(aninfo{2}, animal, 'ca1noisekons');
    end
end

%% plot
if plotfigs
    Pp = load_plotting_params({'riptrigall'});
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        hwin = floor(length(lfpstack(ian).data{1}(1,:,1))/2);
        time = [1/Fp.srate*-hwin:1/Fp.srate:0 1/Fp.srate:1/Fp.srate:1/Fp.srate*hwin];
        for t = 1:length(lfpstack(ian).lfptypes)+1
            gridspacemat(t,:,:) = Pp.tscale(t)*permute(lfpstack(ian).ntrodes*ones(1, ...
                length(time)), [3 2 1]);
        end
%         iderip = 0;
        andef = animaldef(animal);
        tetinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');
        reftet_keys = evaluatefilter(tetinfo, 'isequal($area, ''ref'')');
        ca1tet_keys = evaluatefilter(tetinfo, 'isequal($area, ''ca1'')');
        % make a lookup for day epoch irip irip within dayep
        rel_irip = cell2mat(cellfun(@(x) [1:x]', num2cell(lfpstack(ian).numrips_perep),'un', 0));
        iripmat = [[1:length(lfpstack(ian).epoch)]' lfpstack(ian).day ...
            lfpstack(ian).epoch rel_irip];
        
        for irip = 1:length(iripmat(:,1))
            day = lfpstack(ian).day(irip);
            epoch = lfpstack(ian).epoch(irip);
            ripstarttime = lfpstack(ian).ripStartTime(irip);
            reftet = reftet_keys(ismember(reftet_keys(:,[1 2]),[day epoch], 'rows'),3);
            ca1tets = ca1tet_keys(ismember(ca1tet_keys(:,[1 2]), [day epoch], 'rows'),3); 
            % std thresh
            ripind = knnsearch(rips{ian}{day}{epoch}{1}.starttime,ripstarttime);
            ripstd = rips{ian}{day}{epoch}{1}.maxthresh(ripind);
            if ripstd < Pp.stdthresh
                continue
            end
            % exclude events within 1 second of a noise event
            noiseind = knnsearch(noise{ian}{day}{epoch}{1}.starttime,ripstarttime);
            if abs(noise{ian}{day}{epoch}{1}.starttime(noiseind) - ripstarttime) < 1
                continue
            end
            
            %% ---- init fig----
            if savefigs && ~pausefigs
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            %% PLOT SPIKES
            
%             ha = tight_subplot(2, length(dstack(ian).lfptypes), [.05 .02],[.1 .1],[.05 .05]);
            subaxis(2,Pp.nrow,1, 'Spacing', Pp.Spacing, ...
                'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
            % need to convert irip (across days) into a dayepoch index to
            % get spike data from the cell array.
            
            irelrip = iripmat(irip, end);
%             try % error and reset counter once current day eps are depleted
            irip_muspikes = squeeze(spikestack(ian).muspikes{day}{epoch}(irelrip,:,:))';
%             catch
%                 iderip = 1;
%                 irip_muspikes = squeeze(spikestack(ian).muspikes{day}{epoch}(iderip,:,:))';
%             end
            fprintf('%s %d %d allrip%d derip%d \n', animal, day, epoch, irip, irelrip);
            iripmuntrodes = spikestack(ian).mu_de_keys{day}{epoch}(:,3);
%             irip_muspikes(irip_muspikes>0) = 1;
%             f1 = imagesc(time, iripmuntrodes, irip_muspikes);
%             colormap(1-gray);
            
            [sx, sy] = find(irip_muspikes');
            ca1y = find(ismember(iripmuntrodes, ca1tets));
            sz = zeros(length(sy), 3);
            sz(ismember(sy, ca1y), 3) = 1; % make ca1 blue
            refy = find(ismember(iripmuntrodes, reftet));
            sz(sy==refy, 2) = 1; % make ref green
            
            f1 = scatter(sx/1000-1.001, sy, 10, sz,'+'); % could ad 'z' data for coloring
%             f1.MarkerFaceAlpha = 0.1;
            f1.MarkerEdgeAlpha = 0.05;
%             f1.MarkerEdgeColor = [0 0 0];
            axis tight
            xlim([-1 1])
            yticks(iripmuntrodes);
            yticklabels({num2str(iripmuntrodes)});
            ylabel('ntrode','FontSize',8,'FontWeight','bold', 'FontName','Arial')
            title('mu', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')

            %% plot rip patches in window
            winstarttime = lfpstack(ian).ripStartTime(irip)-(hwin*1/Fp.srate);
            winendtime = lfpstack(ian).ripStartTime(irip)+(hwin*1/Fp.srate);
            ripindsinwin = find(rips{ian}{day}{epoch}{1}.starttime>winstarttime & ...
                rips{ian}{day}{epoch}{1}.endtime<winendtime);
            ripsinwinTimes = [rips{ian}{day}{epoch}{1}.starttime(ripindsinwin) ...
                rips{ian}{day}{epoch}{1}.endtime(ripindsinwin)];
            xs = ripsinwinTimes(:,1)-lfpstack(ian).ripStartTime(irip);
            xe = ripsinwinTimes(:,2)-lfpstack(ian).ripStartTime(irip);
%             xe = xe*1000;
            axis tight
            ye = ylim;
            ye = ye(2);
            %                 ye = max(max(iripdata));
            patch([xs'; xe'; xe'; xs'], repmat([0 0 ye ye]', 1, length(xe)),'y', ...
                'FaceAlpha',.15, 'edgecolor','none');
          %% mu FR
            subaxis(2,Pp.nrow,2, 'Spacing', Pp.Spacing, ...
                'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
            
            irip_muinstaFR = squeeze(spikestack(ian).muinstantFR{day}{epoch}(irelrip,:,:))';
            imagesc(time, iripmuntrodes, 10*log10(flipud(irip_muinstaFR)));
            yticks(iripmuntrodes);
            yticklabels({num2str(flipud(iripmuntrodes))});
            colormap(jet)
%             xlabel('time')
            title('mu firing rate', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
            patch([xs'; xe'; xe'; xs'], repmat([0 0 ye ye]', 1, length(xe)),'k', ...
                'FaceAlpha',0, 'EdgeAlpha', .5, 'edgecolor','k', 'LineWidth', .05);

            % plot sing unit, .. one row for each single unit, labeled or
            % colored by area, and labeled with which ntrode each su came
            % from.
            
            %% plot LFP
            for t = 1:length(lfpstack(ian).lfptypes)
%                 axes(ha(t));
                hold off
                subaxis(2,Pp.nrow,t+2, 'Spacing', Pp.Spacing, ...
                    'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                    Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                    'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                hold on
                if t > 2 % zscore only the bandpass filtered traces
                    iripdata = squeeze(gridspacemat(t,:,:) + ...
                        zscore(double(lfpstack(ian).data{t}(irip, :, :)), [], 2))';
                else
                    try
                        iripdata = squeeze(gridspacemat(t,:,:) + ...
                            lfpstack(ian).data{t}(irip, :, :))';
                    catch
                        iripdata = squeeze(int16(gridspacemat(t,:,:)) + ...
                            lfpstack(ian).data{t}(irip, :, :))';
                    end
                end
                p = plot(time, iripdata);

                ca1y = find(ismember(lfpstack(ian).ntrodes, ca1tets));
                sz = zeros(length(lfpstack(ian).ntrodes), 3);
                sz(ismember(lfpstack(ian).ntrodes, ca1y), 3) = 1; % make ca1 blue
                refy = find(ismember(lfpstack(ian).ntrodes, reftet));
                sz(lfpstack(ian).ntrodes==refy, 2) = 1; % make ref green
                set(p, {'color'}, num2cell(sz, 2));
                
                yticks([squeeze(gridspacemat(t,1,:))]);
                yticklabels({num2str(lfpstack(ian).ntrodes)});
                title(sprintf('%s %s Hz',Fp.LFPtypes{t}, Fp.LFPrangesHz{t}), ...
                    'FontSize',8,'FontWeight','bold', 'FontName','Arial');
                axis tight
                ye = ylim;
                ye = ye(2);
                patch([xs'; xe'; xe'; xs'], repmat([0 0 ye ye]', 1, length(xe)),'y', ...
                    'FaceAlpha',.15, 'edgecolor','none');
            end
            ye = ylim;
            ye = ye(2);
            patch([xs'; xe'; xe'; xs'], repmat([0 0 ye ye]', 1, length(xe)),'y', ...
                'FaceAlpha',.15, 'edgecolor','none');
            for i = 1:length(ripindsinwin)
                text(xe(i), 1, ...
                    sprintf('%.02f',rips{ian}{day}{epoch}{1}.maxthresh(ripindsinwin(i))), ...
                    'FontName', 'Arial', 'FontSize',6);
            end
           %% single unit
            subaxis(2,Pp.nrow, [9 10], 'Spacing', Pp.Spacing, ...
                'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
            
            irip_suspikes = squeeze(spikestack(ian).suspikes{day}{epoch}(irelrip,:,:))';
            iripsu_ntcl = spikestack(ian).sukeys{day}{epoch}(:,[3 4]);
%             irip_muspikes(irip_muspikes>0) = 1;
%             f1 = imagesc(time, iripmuntrodes, irip_muspikes);
%             colormap(1-gray);
            [sx, sy] = find(irip_suspikes');
            ca1y = find(ismember(iripsu_ntcl(:,1), ca1tets));
            sz = zeros(length(sy), 3);
            sz(ismember(sy, ca1y), 3) = 1;
            f1 = scatter(sx/1000-1.001, sy,10,sz,'+'); % could ad 'z' data for coloring
%             f1.MarkerFaceAlpha = 0.1;
            f1.MarkerEdgeAlpha = .4;
%             f1.MarkerEdgeColor = [0 0 0];
            yticks(1:length(iripsu_ntcl(:,1)));
%             yticklabels(, 'FontSize', 6);
            set(gca,'YTickLabel',{num2str(iripsu_ntcl)},'fontsize',6)
%             set(gca,'YTickLabelMode','auto')
            ylabel('ntrode - cluster','FontSize',8,'FontWeight','bold', 'FontName','Arial')
            xticks([-1 -.5 0 .5 1]);
            xlabel('time','FontSize',8,'FontWeight','bold', 'FontName','Arial');
            title('su', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
            axis tight
            xlim([-1 1])
            ye = ylim;
            ye = ye(2);
            patch([xs'; xe'; xe'; xs'], repmat([0 0 ye ye]', 1, length(xe)),'y', ...
                'FaceAlpha',.15, 'edgecolor','none');
            
            %% 2d position
%             axes(ha(t+3));
            subaxis(2,Pp.nrow,[11 12], 'Spacing', Pp.Spacing, ...
                'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
            
            posidx = knnsearch(pos{ian}{day}{epoch}.data(:,1), ...
                lfpstack(ian).ripStartTime(irip));
            posirip = pos{ian}{day}{epoch}.data(posidx,:);
            xpos = pos{ian}{day}{epoch}.data(:,6);
            ypos = pos{ian}{day}{epoch}.data(:,7);
            p = scatter(xpos, ypos, 1, 'filled');
            p.MarkerFaceColor = [0 0 0];
            p.MarkerFaceAlpha = .1;
            hold on
            
            s = posidx - 300;
            if s<1; s=1; end
            e = posidx + 300;
            if e > length(xpos(:,1)); e = length(xpos(:,1)); end
            
            ixpos = xpos(s:e,1);
            iypos = ypos(s:e,1);
            jt = cool(length(iypos));
            be = scatter(ixpos, iypos, 10, jt);
%             af = scatter(pos{ian}{day}{epoch}.data(posidx:e,6), ...
%                 pos{ian}{day}{epoch}.data(posidx:e,7), 5, 'filled');
%             be.MarkerFaceColor = [0 .5 1];
%             be.MarkerFaceAlpha = 0;
%             af.MarkerFaceColor = [.5 0 1];
%             af.MarkerFaceAlpha = .5;
            plot(posirip(6), posirip(7),'.r', 'MarkerSize', 30)
            title('position', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
            xlabel('cm');
%             ylabel('cm');
            axis tight
            xl = xlim;
            yl = ylim;
            if xl(2) > 250 || xl(1) < 0
                xlim([50 250])
            end
            if yl(1) > 250 || yl(1) < 0
                ylim([20 150])
            end
            linposidx = knnsearch(linpos{ian}{day}{epoch}.statematrix.time(:,1), lfpstack(ian).ripStartTime(irip));
            wellexitenter = linpos{ian}{day}{epoch}.statematrix.wellExitEnter(linposidx,:);
%             linpos{ian}{day}{epoch}.statematrix.segmentIndex(linposidx,1)
%             sc = scatter(linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter,1) ...
%                 ,linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter,2), '.');
            wsx = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(1),1);
            wsy = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(1),2);
            wex = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(2),1);
            wey = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(2),2);
            text(wsx, wsy,'\wedge','FontSize',20,'FontWeight','bold', ...
                'Color', [0 .6 0], 'horizontalAlignment', 'center')
            text(wex, wey,'\vee','Fontsize',20,'FontWeight','bold',...
                'Color', [.6 .3 .6], 'horizontalAlignment', 'center')
%             linpos{ian}{day}{epoch}.statematrix.segmentIndex(linposidx,:);
            
            hold off
            
            %% Linear position
            subaxis(2,Pp.nrow,[13 16], 'Spacing', Pp.Spacing, ...
                'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
%             postime = pos{ian}{day}{epoch}.data(:,1);
%             scatter(postime, xpos, abs(zscore(ypos).^10),xpos) % this is too fancy
%             colormap(hsv);
            
            postime = linpos{ian}{day}{epoch}.statematrix.time;
            traj = linpos{ian}{day}{epoch}.statematrix.traj;
            segmentnum = linpos{ian}{day}{epoch}.statematrix.segmentIndex;
%             % plot animal's travel lindist from center well and colorize by segment
            lindist = linpos{ian}{day}{epoch}.statematrix.linearDistanceToWells(:,1);
            % plot linpos
            plot(postime,lindist,'linewidth',1,'color',[.8 .8 .8]);
            hold on
            
            %         clr1 = [.5 .3 .8; .2 .1 .1; .4 .1 .6; .2 .9 .6];
            clr1 = [.6 .5 .2; .8 .5 .8; .6 0 .6; 0 0.4980 0.3255; 0 0.8 0.5216]; %.2 .9 .6];
            %             %plot ripple times from all areas
            for iseg = 1:5
                inds = [];
                inds = double(segmentnum == iseg);
                inds(find(inds==0))=nan;
                plot(postime.*inds,lindist.*inds,'-','linewidth',2,'color',clr1(iseg,:));
                hold on;
            end
           %% performance patches over linpos
            trialIO = behave{ian}.BehaveState.statechanges{day}{epoch}.statechangeseq;
            outBcol = 9; lastTimecol= 2; currTimecol=3; corrcol=7; timeportoutcol=10;
            outBCorr = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==1), ...
                [lastTimecol currTimecol corrcol timeportoutcol]);
            outBMist = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==0), ...
                [lastTimecol currTimecol corrcol timeportoutcol]);
            inBCorr = trialIO((trialIO(:,outBcol)==0 & trialIO(:,corrcol)==1), ...
                [lastTimecol currTimecol corrcol timeportoutcol]);
            inBMist = trialIO((trialIO(:,outBcol)==0 & trialIO(:,corrcol)==0), ...
                [lastTimecol currTimecol corrcol timeportoutcol]);
            clr2 = ([.5 .5 1; .5 .7 .7; .3 .8 .3; .5 .5 1; .5 .7 .7; .3 .8 .3;]);
            yl = ylim;
            xl = xlim;
            for itrial = 1:length(outBCorr(:,1)) %patch correct outbound trials
                patch([outBCorr(itrial,1) outBCorr(itrial,2) outBCorr(itrial,2) ...
                    outBCorr(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'y', 'edgecolor', ...
                    'none', 'FaceAlpha',.2);
            end
            for itrial = 1:length(outBMist(:,1)) %patch incorrect outbound trials
                patch([outBMist(itrial,1) outBMist(itrial,2) outBMist(itrial,2) ...
                    outBMist(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'k', 'edgecolor', ...
                    'none', 'FaceAlpha',.1);
            end
            for itrial = 1:length(outBCorr(:,1)) %is a portOUT
                plot(outBCorr(itrial,4), yl(2)-1, 'd', 'linewidth', 1, 'MarkerEdgeColor', ...
                    'k','MarkerFaceColor', 'y'); hold on;
            end
            % inbound patches
            for itrial = 1:length(inBCorr(:,1)) %patch correct outbound trials
                patch([inBCorr(itrial,1) inBCorr(itrial,2) inBCorr(itrial,2) ...
                    inBCorr(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'y', 'edgecolor', ...
                    'none', 'FaceAlpha',.2);
            end
            for itrial = 1:length(inBMist(:,1)) %patch incorrect outbound trials
                if inBMist(itrial,1) == 0
                    continue
                end
                patch([inBMist(itrial,1) inBMist(itrial,2) inBMist(itrial,2) ...
                    inBMist(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'k', 'edgecolor', ...
                    'none', 'FaceAlpha',.1);
            end
            for itrial = 1:length(inBCorr(:,1)) %is a portOUT
                plot(inBCorr(itrial,4), yl(1)-2, 'd', 'linewidth', 1, 'MarkerEdgeColor', ...
                    'k','MarkerFaceColor', 'y'); hold on;
            end
            
             % plot all rips in ep over linpos
            de_ripidx = find(ismember([lfpstack(ian).day lfpstack(ian).epoch], ...
                [day epoch], 'rows'));
            de_ripidxlinpos = knnsearch(postime, lfpstack(ian).ripStartTime(de_ripidx));
            f = scatter(lfpstack(ian).ripStartTime(de_ripidx), lindist(de_ripidxlinpos), ...
                10, 'filled');
            f.MarkerFaceColor = 'k';
            % highlight current rip
            hold on
            linposripidx = knnsearch(postime, lfpstack(ian).ripStartTime(irip));
            plot(lfpstack(ian).ripStartTime(irip), lindist(linposripidx), ...
                '.r', 'MarkerSize', 30)
            
            title('linpos + performance', 'FontSize',8,'FontWeight','bold', 'FontName', ...
                'Arial')
            xlabel('time s');
%             ylabel('xpos cm');
            axis tight
%             %% all performance
%             subaxis(2,7,[13 14], 'Spacing', Pp.Spacing, ...
%                 'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
%                 Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
%                 'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
%             
%             trialIOfields = BehaveState.statechanges{day}{ep}.fields;
%             corrcol = find(cell2mat(cellfun(@(x) strcmp(x,'correct'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
%             timeportoutcol = find(cell2mat(cellfun(@(x) strcmp(x,'timeportout'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
%             outBcol = find(cell2mat(cellfun(@(x) strcmp(x,'outbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
%             inBcol = find(cell2mat(cellfun(@(x) strcmp(x,'inbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
%             lastTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'lasttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
%             currTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'currenttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            
%             
%             axis tight
%             set(gca,'TickDir','out')
%             a = gca;
%             % set box property to off and remove background color
%             set(a,'box','off','color','none')
% 
% %             ssptimes = timesbehave{ian}.BehaveState.statespace.allepsmat(:,3);
% %             ssplasttimes = timesbehave{ian}.BehaveState.statespace.allepsmat(:,3);
% %             
% %             behave{ian}.BehaveState.sequences.
% %             
% %             
%             
%             % plot instead the trajectory since the last well visit, and
%             % the next well visit, and whether they were rewarded or not, and which ones.
%             
%             % behave state context (where on the total
%             % learning curve the current rip is, if he is in an incorrect/correct,
%             % outbound/inbound state.
            
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %d %d %s %s rip%d start%dms', animal, day, epoch, ...
                Fp.paths.filenamesave(5:end), Fp.epochEnvironment, irip, ...
                round(1000*lfpstack(ian).ripStartTime(irip)));
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 12);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(Fp.paths.figdirectory, 'riptrig_all', sprtit)
                close all
            end
        end
    end
end