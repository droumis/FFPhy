
        

            figure; hold on
            
            % Spike raster plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,3,[1 4 7]);
            
            % plot 2 SD SWRs that occur in the widnow %%%%%%%%%
            consvec = output1{day}{ep}.cons;
            c = lookup(ev_starttime - 0.5,consvectimes);
            d = lookup(ev_starttime + 0.5,consvectimes);
            if 0
                plot(1000 * (consvectimes_rip2(c:d) - rip_starttime), consvec_rip2(c:d) + 0.5,'-','Color',[1 .7 .7],'linewidth',8);      hold on
            else
                % identify which 2 SD ripples to plot
                    % first find 2 SD ripples again (brief calculation)
                    ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep ],'ripplescons',1,...
                        'consensus_numtets',3,'minthresh',2,...
                        'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
                    consvec_rip2 = ripout{day}{ep}.cons;
                    consvectimes_rip2 = ripout{day}{ep}.time;
                    periodtimes_rip2 = vec2list(consvec_rip2,consvectimes_rip2);
                plotwinvec = logical(list2vec([ev_starttime-0.5   ev_starttime+0.5],ripout{day}{ep}.time))';
                consvec_rip2_toplot = consvec_rip2 & plotwinvec;
                rip_toplot = 1000 * (vec2list(consvec_rip2_toplot,ripout{day}{ep}.time) - ev_starttime);
                numrip_toplot = size(rip_toplot,1);
                % plot all 2 SD ripples in the raster window
                for rp = 1:numrip_toplot
                    patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)],...
                        [1 1 (maxtet + 4) (maxtet + 4)],...
                        [1 .9 .9],'edgecolor','none'); hold on
                end
            end
            % plot EVENT as a thick red line (i.e. ripple or wave gamma) %%%%%%%%%%%
            hold on
            plot([ev_startms ev_endms] - ev_startms,[.5 .5],'Color','r','linewidth',4)
            if 0
                set(gca,'xtick',[-300:100:300]);
                set(gca,'xticklabel',{'-300','','','0','','','+300'});
                xlim([-300 300]);
            else
                set(gca,'xtick',[-500:100:500]);
                set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
                xlim([-500 500]);
            end
            
            % plot SPIKES raster %%
            for tet = 1:maxtet
                spikebins = find(spike_mats(tet,:,ww) > 0);
                numlines = length(spikebins);
                for ll = 1:numlines
                    plot([spikebins(ll) spikebins(ll)] - extratime,[tet-0.5 tet+0.5],...
                        'linestyle','-','Color',[0 0 0],'LineWidth',1)
                end
            end
            set(gca,'ytick',1:maxtet)
            %set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
            xlabel('Time (ms)') ;   ylabel('Tetrode');
            title([cons_name1 ' # : ',num2str(ww)],'fontweight','bold','fontsize',14);
            set(gca,'tickdir','out');
            
            box off
            
            % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
            a = lookup(ev_starttime - 0.5,riptrace_timevec) ;
            b = lookup(ev_starttime + 0.5,riptrace_timevec) ;
            plot(1000 * (riptrace_timevec(a:b) - ev_starttime), -riptrace(a:b)/2 + maxtet + 4,'r-','linewidth',2); hold on
            a = lookup(ev_starttime - 0.5,wgtrace_timevec) ;
            b = lookup(ev_starttime + 0.5,wgtrace_timevec) ;
            plot(1000 * (wgtrace_timevec(a:b) - ev_starttime), 2 * -wgtrace(a:b)/2 + maxtet + 4,'-','Color',[0 .5 1],'linewidth',2);
            
            set(gca,'ydir','reverse')
            set(gca,'fontsize',12)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Decoded posterior image plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for n = 1:num_encodeeps 
                if n == 1
                    subplot(3,3,[2 5 8]);
                elseif n == 2
                    subplot(3,3,[3 6 9]);
                end
                linposbins = (1:length(xbins_cut{n})) * xdel ;
                imagesc(xvecms,linposbins,posteriors{n}(:,:,end),[0 .1]);
                
                %title('postx, clusterless','FontSize',12);
                ylabel('Linearized position','FontSize',12);
                %set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
                %colormap(flipud(hot(256)));
                colormap(hot);
                caxis([0 0.1]);
                ylim([-.01 max(linposbins)])
                set(gca,'tickdir','out');
                
                set(gca,'xtick',[-500:100:500]);
                set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
                xlim([-500 500])
                set(gca,'fontsize',12)
                xlabel('Time (ms)')
                
                hold on
                % in local environment, plot current position (thick grey line)
                if n == 1
                    winpos = armdists_cat{n}( ev_posstartind : ev_posendind );  % armdists_cut (linear) positions over the course of the ripple
                    winpos_timevec = 1000 * ( postimevec{1}(ev_posstartind:ev_posendind) - ev_starttime );  % 1-ms indices over the course of the ripple
                    plot(winpos_timevec,winpos,'-','linewidth',6,'Color',[.8 .8 .8])
                end
                
                % plot line indicating start and end of ripple
                %plot([ripms_start ripms_end] - ripms_start,[50 50],'Color','r','linewidth',4)
                plot([0 0],[0 max(linposbins)],'Color','r','linewidth',1)
                plot([ev_endms - ev_startms ev_endms - ev_startms],...
                    [0 max(linposbins)],'Color','r','linewidth',1)
                
                % plot lines demarcating positions of arms
                % between Center and Right
                plot([-500 500],[centerarmmax(n) centerarmmax(n)],'--','Color','w','linewidth',2)
                % between Right and Left
                plot([-500 500],[centerarmmax(n)+rightarmmax(n) ...
                    centerarmmax(n)+rightarmmax(n)],'--','Color','w','linewidth',2)
            end
            


