%{

XP-triggered average of z scored ripple power trace
as suggested by david
Meeting Notes: https://docs.google.com/document/d/1wTUtn20BHk7W64m5979wIZ19XYAZmgIi7sacT9YJoQc/edit#bookmark=id.fi6y0ntpmyqv

copy how i got z rip pwr from lickswrExamples_20190926
%}


pconf = paramconfig;
create_filter = 1;
run_ff = 0;
load_ffdata = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'pdf', 'png'};

plot_XPtrigAvgRip_pAn_pDay = 0;
plot_XPtrigAvgRip_pAn = 1;
plot_XPtrigAvgRip = 0;
%% FF Data
Fp = [];
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ4'};
Fp.filtfunction = 'dfa_XPtrigAvgRip'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'XPtrigAvgRip';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    'proximalWell', '<4cm/s', Fp.Label, Fp.filtfunction};

Fp = load_filter_params(Fp);

if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,...
        'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F),...
        'un', 1);
    F = runfilter(F);
    
%     save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
%         'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_' Fp.Label]);
end

%% Plot per day for each animal. then per animal. then all
if plotfigs
    if plot_XPtrigAvgRip_pAn_pDay
        figname = sprintf('%s-pAn-pDay','xcNormAC');
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            days = cell2mat([{F(a).output{1}.index}']);
            ndays = size(days,1);
            ncols = 1;
            for d = 1:ndays
                day = days(d,1);
                idata = F(a).output{1}(d);
                % mean trace per an per day
                sf1 = subaxis(ndays, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
%                 ripPowPredCorrCoef = nanmean(idata.ripPowPredCorrCoef);

                
%                 xcnorm = idata.xcnorm;
%                 semd = std(xcnorm,[],1)/sqrt(size(xcnorm,1));
                t = -.5:1/1500:.5;
                xcnorm = idata.xcnorm;
%                 fill([t'; flipud(t')],[xcnorm-semd;flipud(dm+semd)], 'b',...
%                     'linestyle', 'none', 'facealpha', .2);
%                 hold on
                plot(t, xcnorm, 'color', [0 0 1 1], 'linewidth', 2);
                ylim([-.5 .5])
                line([0 0], ylim, 'linestyle', '--', 'color', 'k')
                
%                 text(-.45,.45,sprintf('ripPowPredCorrCoef %f', ripPowPredCorrCoef), 'fontsize', 16)
                stit = sprintf('%s wF%d %s', 'xcNormAC', 0, animal); %Fp.wienerFiltDiv
            
%                 t = idata.time;
%                 st = knnsearch(t', Pp.win(1));
%                 en = knnsearch(t', Pp.win(2));
%                 t = t(st:en);
% %                 m = idata.mean_rippwr_XPtrig;
%             
%                 m = idata.mean_rippwr_XPtrig(:,st:en);
%                 s = idata.sem_rippwr_XPtrig(:,st:en);
% %                 mall = nanmean(nanmean(rippwr_XPtrig));
% %                 rippwr_XPtrig_mcenter = rippwr_XPtrig - mall;
% %                 m = nanmean(rippwr_XPtrig_mcenter);
% %                 s = nanstd(rippwr_XPtrig_mcenter,1)/sqrt(size(rippwr_XPtrig_mcenter,1));
% %                 s = idata.sem_rippwr_XPtrig;
%                 
% 
% %                 m = m(st:en);
% %                 s = s(st:en);
%                 
%                 fill([t'; flipud(t')],[m'-s';flipud(m'+s')], 'k',...
%                     'linestyle', 'none', 'facealpha', .4);
%                 hold on
%                 plot(t, m, 'color', [0 0 0 1], 'linewidth', 1);
%                 title(sprintf('day %d', d));
%                 ylabel('rippwr')
%                 xticks(Pp.win(1):.5:Pp.win(2))
%                 xlabel('time from XP');
%                 axis tight
%                 line(xlim, [0 0], 'linestyle', '-', ...
%                     'color', [.5 .5 1 .5], 'linewidth', 2)
%                 ylim([-.75 .75])

                
            end
            xlabel('time from BB')
            % super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname],...
                    stit, 'savefigas', savefigas);
            end
        end
    end
    if plot_XPtrigAvgRip_pAn
        figname = sprintf('%s-pAn',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            %%
            if 0
                shd = horzcat(F(a).output{1}.circperm_xcnorm);
                shd = shd(:)';
                shd = cell2mat(shd);
                m = nanmean(shd,2);
                s = nanstd(shd,[],2);

                fill([t'; flipud(t')],[m-s;flipud(m+s)], 'k',...
                    'linestyle', 'none', 'facealpha', .1);
                hold on;
                fill([t'; flipud(t')],[m-s*2;flipud(m+s*2)], 'k',...
                    'linestyle', 'none', 'facealpha', .2);
            end
            ripPowPredCorrCoef = nanmean(vertcat(F(a).output{1}.ripPowPredCorrCoef));
            
            
            d = horzcat(F(a).output{1}.xcnorm);
            semd = std(d,[],2)/sqrt(size(d,2));
            t = -.5:1/1500:.5;
            dm = nanmean(d,2);
            fill([t'; flipud(t')],[dm-semd;flipud(dm+semd)], 'b',...
                'linestyle', 'none', 'facealpha', .2);
            hold on
            plot(t, dm, 'color', [0 0 1 1], 'linewidth', 2);
            ylim([-.5 .5])
            line([0 0], ylim, 'linestyle', '--', 'color', 'k')
            xlabel('time from BB')
            text(-.45,.45,sprintf('ripPowPredCorrCoef %f', ripPowPredCorrCoef), 'fontsize', 16)
            stit = sprintf('%s wF%d %s', 'xcNormAC', 0, animal); %Fp.wienerFiltDiv
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' 'xcNormAC'],...
                    stit, 'savefigas', savefigas);
            end
        end
    end            
            %%
%             t = F(a).output{1}(1).time;
%             st = knnsearch(t', Pp.win(1));
%             en = knnsearch(t', Pp.win(2));
%             t = t(st:en);
%             fs = 1500;
            % plot xcorr
%             subaxis(3, 2, 1, Pp.posparams{:}); 
%             plot(d)
%             d = d(1:end-1,:);
%             for x = 1:size(d,1)
%                 plot(d(x,:)); hold on
%                 fprintf('%d\n', x)
%             end
%             legend; 
%             m_ifft_xc = nanmean(zscore(d,[],2),1);
%             m_ifft_xc = nanmean(d);
%             plot(t, d, 'linewidth', 2); 
% %             line([0 0], ylim, 'color', 'k')
%             axis tight;
%             xlabel('time')
% %             imagesc(real(ifft_xc_acdiv_XPtrig));
%             
%             
% %             xc = nanmean(vertcat(F(a).output{1}.mean_rippwr_XPtrig),1);
% %             xc = xc(:,st:en);
% %             s = nanmean(vertcat(F(a).output{1}.sem_rippwr_XPtrig),1);
% %             s = s(:,st:en);
% %           
%             plot(t, xc, 'color', [0 0 0 1], 'linewidth', 1);
%             hold on;
%             fill([t'; flipud(t')],[xc'-s';flipud(xc'+s')], 'k',...
%                 'linestyle', 'none', 'facealpha', .2);            
%             
%             ylabel('rippwr')
%             xticks(Pp.win(1):.5:Pp.win(2))
%             title('xcorr');
%             axis tight;
%             hold off
%             ylim([-.5 .5])
%             
%             % plot xcorr spectrum
%             subaxis(3, 2, 2, Pp.posparams{:});
%             [xcPxx,Freq] = periodogram(xc,rectwin(length(xc)),length(xc),fs);
%             plot(Freq,10*log10(xcPxx), 'linewidth', 2)
%             grid on;
%             xlabel('Hz'); ylabel('dB/Hz');
%             title('fft(xcorr)');
%             xlim([0 20])
%             ylim([-100 0])
%             
%             % plot XP acorr
%             subaxis(3, 2, 3, Pp.posparams{:});
%             XPacorr = horzcat(F(a).output{1}.XPacorr)';
%             XPacorrm = nanmean(nanmean(XPacorr));
%             m = nanmean(XPacorr-XPacorrm,1);
%             m = m(:,st:en);
%             
% %             m(lookup(0,t)) = 0;
%             plot(t, m, 'color', [0 0 0 1], 'linewidth', 2);
%             title('acorr')
%             axis tight
%             ylim([-.008 .008])
% 
%             % plot acorr spectrum
%             subaxis(3, 2, 4, Pp.posparams{:});
%             [acPxx,Freq] = periodogram(m,rectwin(length(m)),length(m),fs);
%             plot(Freq,10*log10(acPxx), 'linewidth', 2)
%             grid on;
%             xlabel('Hz'); ylabel('dB/Hz');
%             title('fft(acorr)');
%             xlim([0 20])
%             ylim([-100 0])
%             
%             
%             % plot normed xcorr
%             subaxis(3, 2, 5, Pp.posparams{:});     
%             xcfft = fft(xc);
%             acfft = fft(m);
%             f = xcfft ./ acfft;
%             s = real(ifft(f));
%             plot(t, s, 'color', [0 0 0 1], 'linewidth', 2);
%             title('ifft(fft(xc) / fft(ac))')       
%             yl = ylim;
%             line([0 0], yl, 'linestyle', '--', 'linewidth', 1, 'color', 'm')
%             
%             % plot normed xcorr zoomed in
%             subaxis(3, 2, 6, Pp.posparams{:});     
%             plot(t, s, 'color', [0 0 0 1], 'linewidth', 2);
%             axis tight
%             xlim([-.5 .5])
%             title('ifft(fft(xc) / fft(ac)) zoomX')       
%             yl = ylim;
%             line([0 0], yl, 'linestyle', '--', 'linewidth', 1, 'color', 'm')
%             
% %             % plot norm xcorr spectrum
% %             subaxis(3, 2, 6, Pp.posparams{:});
% %             [xcnPxx,Freq] = periodogram(s,rectwin(length(m)),length(m),fs);
% %             plot(Freq,10*log10(xcnPxx), 'linewidth', 2)
% %             grid on;
% %             xlabel('Hz'); ylabel('dB/Hz');
% %             title('fft(xcorr');
% %             xlim([0 20])

            % super

    if plot_XPtrigAvgRip
        figname = sprintf('%s',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        
        % process
        rippwr_XPtrig = [];
        XPpostSWR_offset = [];
        for a = 1:length(F)
            rippwr_XPtrig = [rippwr_XPtrig; vertcat(F(a).output{1}.rippwr_XPtrig)];
        end
        m = nanmean(rippwr_XPtrig,1);
        s = nanstd(rippwr_XPtrig,1)/sqrt(size(rippwr_XPtrig,1));
        t = idata.time;
        st = knnsearch(t', Pp.win(1));
        en = knnsearch(t', Pp.win(2));
        t = t(st:en);
        m = m(st:en);
        s = s(st:en);
        % plot
        ifig = init_plot(showfigs, Pp.position);
        plot(t, m, 'color', [0 0 0 1], 'linewidth', 1);
        hold on;
        fill([t'; flipud(t')],[m'-s';flipud(m'+s')], 'k',...
            'linestyle', 'none', 'facealpha', .2);
        ylabel('rippwr')
        xticks(Pp.win(1):.5:Pp.win(2))
        xlabel('time from XP');
        axis tight
        % super
        stit = sprintf('%s', figname);
        setSuperAxTitle(stit);
        if pausefigs
            pause
        end
        if savefigs
            strsave = save_figure([pconf.andef{4} '/' figname],...
                stit, 'savefigas', savefigas);
        end
    end
end