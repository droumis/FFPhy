
% Plots output of clusterless decode (Pos)

clear all

plotter         = 1;
plot_placemap   = 1;
place_map_only  = 0;
    save_plot   = 0;
    
    EXAMPLE                         = 12;              %    
        clusterless_examples;   % <<< select examples here
        EXTYPE_toplot               = 1;
        FIGURESIZE                  = [.1 .06 .37 .6];
        DECODE_RUN                  = 2;     % 0: Uniform, 1: Randwalk, 2: Randwalk_new
        PLOT_OUTER_INTEGRATION      = 0;
        MAXVAL                      = 0.03;
        GREEN_POS_CHOICE            = 2;
        
        
    if place_map_only
        plotter = 0;
        analyze = 0;
        plot_placemap = 0;
    end
    
    omit_mainplot = 0;
    
    if plotter 
        
        
        MOVINGPERIODS = 2 ;  % 1: vel4, 2: nonimmobile05, 3: nonimmobile1
        
        % Study figures %%%
        plot_Xcorr_LR = 0;
        plot_LRprop_hist = 0;
        plot_LRprop_phasehist = 0;
        plot_altspeed_hist = 1;
        plot_pvalue_scatter = 1;
        plot_altspeed_counts = 1;
        
        % LFP plot %%%%
        LFPMAX = 0.5;  % mV
        plot_LFP = 0;
        plot_LOWFREQ_LFP = 0;      
        plot_DELTA_LFP = 0;  
        
    end
        % place plot %%%%%%%%%%%%%%%%%%%%%%%%%   
        DARKGREY  = [.72 .72 .72];
        ANGVEL_SCALE = 2;
        NUMSAMP_HEADPOS = 2 ;    
        
    % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotter
        
        if DECODE_RUN == 0
            decodedir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Uniform';
            thetabindir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Uniform/Thetabins_LR';
        elseif DECODE_RUN == 1
            decodedir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Randwalk';
            thetabindir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Randwalk/Thetabins_LR';
        elseif DECODE_RUN == 2
            decodedir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Randwalk_new/';            
            thetabindir = '';
        end
                
        disp(['DECODE_RUN: ' num2str(DECODE_RUN)])
        
        cd(decodedir)
        
        % loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
            'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
        animpref = animal_toplot(1:3);
        animalinfo = animaldef(animal_toplot);
        daydir = getdaydir(animal_toplot);
        an = find(strcmp(animal_toplot,animals_order));
        
        d = dayep(1);
        ep = dayep(2);

        % load Decode file
        if ~exist('P','var') || isempty(P) || ~strcmp(P.animalname,animal_toplot) || ...
               ~all(P.dayep == dayep) || P.EXTYPE ~= EXTYPE_toplot
           
           % Decode file %%%%%%%%%%%%%%
           filestr = sprintf('%s_%d_%d_Pos_%d.mat',animpref,d,ep,EXTYPE_toplot);
           disp(sprintf('loading %s',filestr))
           load(filestr,'P');
            
           % Positional %%%%%%%%%%%%%%%%%%%
           pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',P.dayep(1));
            [veltrace,hdir,xpos,ypos] = posparser(pos{d}{ep}.data);
           linpos = loaddatastruct(animalinfo{2},animalinfo{3},'linpos',P.dayep(1));
           armpos = loaddatastruct(animalinfo{2},animalinfo{3},'armpos',P.dayep(1));
           lindist = linpos{d}{ep}.statematrix.lindist;
                CPbuff = choicepointer(linpos{d}{ep});
                CPprevec = lindist < CPbuff;
                CPvec = (lindist < CPbuff );   % 10 cm zone  (in green)  %
                Lvec =  ismember(linpos{d}{ep}.statematrix.segmentIndex,[4 5]);
                Rvec =  ismember(linpos{d}{ep}.statematrix.segmentIndex,[2 3]);
                
           postimevec = pos{d}{ep}.data(:,1);    
            epstart = postimevec(1);
            epend = postimevec(end);

            % tets, clust spikes
            spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
            cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
            tetinfo = loaddatastruct(animalinfo{2},animalinfo{3},'tetinfo');
            [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
            maxcell = size(adtc_list,1);
            selected_tets = unique(P.selected_tets);
                numtets = length(selected_tets);           
           
           % LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetabins = loaddatastruct(animalinfo{2},animalinfo{3},'thetabins',d);
            thetatet = thetabins{d}.thetatet;
            if strcmp(animpref,'Dav') && d == 3
                disp('Dav thetatet = 1, d 3 ep 2');
                thetatet = 1;
            end
            theta = loadeegstruct(animalinfo{2},animalinfo{3},'theta',d,ep,thetatet);

            if plotter
                if plot_LFP
                    eeg = loadeegstruct(animalinfo{2},animalinfo{3},'eeg',d,ep,thetatet);
                    eegtrace = eeg{d}{ep}{thetatet}.data;
                    eegtimevec = geteegtimes(eeg{d}{ep}{thetatet});
                end
                if plot_LOWFREQ_LFP
                    lowfreq = loadeegstruct(animalinfo{2},animalinfo{3},'lowfreq',d,ep,thetatet);
                    lowfreqtrace = lowfreq{d}{ep}{thetatet}.data;
                    lowfreqtimevec = geteegtimes(lowfreq{d}{ep}{thetatet});
                end
                if plot_DELTA_LFP
                    delta = loadeegstruct(animalinfo{2},animalinfo{3},'delta',d,ep,thetatet);
                    deltatrace = delta{d}{ep}{thetatet}.data;
                    deltatimevec = geteegtimes(delta{d}{ep}{thetatet});
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           timefilterscript
           if MOVINGPERIODS == 1
               [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{vel4},[d ep]);   
           elseif MOVINGPERIODS == 2
               [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{nonimmobile05},[d ep]);    
           elseif MOVINGPERIODS == 3
               [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{nonimmobile1},[d ep]);   
           end
            movingperiods = dummy{d}{ep};   %
            movingvec = logical(list2vec(movingperiods,postimevec)); 
            
        end
                
        % Theta LFP %
        if strcmp(animpref,'Cor')
            LFPSIGN = 1;
        else
            LFPSIGN = -1;
        end
        thetatrace = double(theta{d}{ep}{thetatet}.data(:,1))/1000;
        thetaphase = double(theta{d}{ep}{thetatet}.data(:,2))/10000;  % hilbert phase
        thetatimevec = geteegtimes(theta{d}{ep}{thetatet});
        
        % Linearized pos %
   
            %posvec_all = P.posvec;
  
            [posvec_all]            = posvecmaker2(lindist,linpos{d}{ep})  ;
            [posvec_all, outersep]  = posvecmaker_compress(posvec_all)     ;
            
       
        % Head pos %
        headang = loaddatastruct(animalinfo{2}, animalinfo{3}, 'headang', d);
            angvel = abs(headang{d}{ep}.angvel)  * 180 / pi;    % convert to degrees

        
        % Identify what to plot (excursions or manual) %%
        
        if ~isempty(plot_start)
            numplots = 1;
%             foundflag = 0;
%             for exx = 1:numexcur
%                 decbins = [ P.Tvecs{exx}(1:(end-1))   P.Tvecs{exx}(2:end) ] ;
%                 if any(logical(isExcluded( plot_start:plot_end , [decbins(1,1) decbins(end,2)] )))
%                     foundflag = 1;
%                     tvec = P.Tvecs{exx};
%                     break
%                 end
%             end
%             if ~foundflag
%                 disp('target period not found')
%                 keyboard
%             end
%             disp(sprintf('%s: day %d ep %d (time: %d - %d) (EXCUR: %d)',P.animalname(1:3),d,ep,numexcur,plot_start,plot_end,EXTYPE_toplot))            
        else
            numplots = length(P.Posts);
            disp(sprintf('%s: day %d ep %d (%d excurs) (EXCUR: %d)',P.animalname(1:3),d,ep,numplots,EXTYPE_toplot))
        end
               
    end
    

% Plotter of posterior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotter

    for PL = 1:numplots
        
        K = figure('units','normalized','outerposition',FIGURESIZE);
        
       % Plot Posterior from each contained excursion %%%%%%%%%%%%%%%

       for xx = 1:length(P.Posts)
           
           tvec = P.Tvecs{xx};
                ex_start = P.Tvecs{xx}(1)   ;
                ex_end   = P.Tvecs{xx}(end) ;
                
            if ~any(isExcluded(ex_start:0.1:ex_end,[plot_start plot_end]))
                continue
            end
           
           decbins = [ tvec(1:(end-1))   tvec(2:end) ];
           a = lookup(plot_start,tvec);
           b = lookup(plot_end,tvec);
           aa = lookup(plot_start,postimevec - postimevec(1));
           bb = lookup(plot_end,postimevec - postimevec(1));
                %posvec2 = posvec(aa:bb);    
                %posvec2(posvec2 <= 0) = nan;
           
           % Process posterior into image
           post        = P.Posts{xx}(:,a:b);
           postfull    = P.Posts{xx}(:,:);
           
           % convert to RGB
           post_rgb = nan(size(post,1),size(post,2),3);  % x, y, (r,g,b)
           yvec = 1:size(post,1);
           Cinds = yvec <= CPbuff;                          % Center arm
           Linds = (yvec > CPbuff) & (yvec <= outersep);    % Left arm
           Rinds = ~Cinds & ~Linds;                         % Right arm
           % C arm
           %maxval = max(post(:));
           maxval = MAXVAL;
           post_rgb(Cinds,:,1) = (maxval - post(Cinds,:))/maxval;
           post_rgb(Cinds,:,2) = (maxval - post(Cinds,:))/maxval;
           post_rgb(Cinds,:,3) = (maxval - post(Cinds,:))/maxval;
           % L arm  (red)
           post_rgb(Linds,:,1) = 1;
           post_rgb(Linds,:,2) = (maxval - post(Linds,:))/maxval;
           post_rgb(Linds,:,3) = (maxval - post(Linds,:))/maxval;
           % R arm  (blue)
           post_rgb(Rinds,:,1) = (maxval - post(Rinds,:))/maxval;
           post_rgb(Rinds,:,2) = (maxval - post(Rinds,:))/maxval;
           post_rgb(Rinds,:,3) = 1;
           
           % saturation
           post_rgb(post_rgb < 0) = 0;
           
           
           
           % Plot image of posterior %%%%%%%%%%%%%%%%%%%%%%%%%
           maxpos = size(post_rgb,1);
           AX(1) = subplot(7,1,[2 3]);
           image(tvec(a:b),1:maxpos,post_rgb); hold on
           colormap jet
           %colormap(flipud(gray))
           set(gca,'ydir','normal')
           set(gca,'xticklabel','')
           xlim([plot_start plot_end])
           
           % Plot LR density trace %%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if PLOT_OUTER_INTEGRATION
               AX(2) = subplot(7,1,6);
               %plot(tvec(a:b),sum(postfull(Cinds,a:b),1),'k-','linewidth',2); hold on
               if 1
                   plot(tvec(a:b),sum(postfull(Linds,a:b),1),'r-','linewidth',2); hold on
                   plot(tvec(a:b),sum(postfull(Rinds,a:b),1),'b-','linewidth',2); hold on
               else
                   smoothL = smoothvect(sum(postfull(Linds,a:b),1),muakern);
                   smoothR = smoothvect(sum(postfull(Rinds,a:b),1),muakern);
                   plot(tvec(a:b),smoothL,'r-','linewidth',2); hold on
                   plot(tvec(a:b),smoothR,'b-','linewidth',2); hold on
               end
               xlim([plot_start plot_end])
               ylim([0 1.1])  ;
               set(gca,'xticklabel','')
               set(gca,'fontsize',14)
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           
       end
       
       % Plot Theta bins on Decodes
       if 0
           cd(thetabindir)
           filestr = sprintf('Thetabins_LR_%d_%s_%d_%d.mat',EXTYPE_toplot,animal_toplot(1:3),d,ep);
           filename = dir(filestr);
           load(filename.name,'tbins_LR');
           tbins = tbins_LR.tbins;
           cd(decodedir)
           
           AX(1) = subplot(7,1,[2 3]);
           Thetabinplotter(tbins-epstart,plot_start,plot_end,[0 maxpos]);
           
           AX(2) = subplot(7,1,6);
           Thetabinplotter(tbins-epstart,plot_start,plot_end,[0 1.1]);
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

    % Plot animal's position on Post image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AX(1) = subplot(7,1,[2 3]);
    % lindist (all)  (green stripe)
    jumpinds = find(   diff([nan posvec_all(:)']) > 30  );  % where goes from center arm to R arm
    posvec_all_2 = posvec_all;
    posvec_all_2(jumpinds) = nan;
    %plot(postimevec(aa:bb) - postimevec(1),posvec_all_2(aa:bb),'k-','linewidth',5);
    plot(postimevec(aa:bb) - postimevec(1),posvec_all_2(aa:bb),'-','Color',[40 200 40]/255,'linewidth',4); hold on
    % plot that single jumped sample!
    for jj = jumpinds(:)'
        linkinds = jj:(jj+1);
        %plot(postimevec(linkinds) - postimevec(1),posvec_all(linkinds),'k-','linewidth',5);
        plot(postimevec(linkinds) - postimevec(1),posvec_all(linkinds),'-','Color',[40 200 40]/255,'linewidth',4); hold on
    end

    if 0
        % posvec (encoding positions, taken from P file) (stripe)
        % identify jumps and put a single nan for sake of plotting
        jumpinds = find(diff([nan posvec(:)']) > 30 );  % where goes from center arm to R arm
        posvec2 = posvec;
        posvec2(jumpinds) = nan;
        plot(postimevec(aa:bb) - postimevec(1),posvec2,'k-','linewidth',5); hold on
        plot(postimevec(aa:bb) - postimevec(1),posvec2,'-','Color',[40 200 40]/255,'linewidth',2); hold on
        % plot that single jumped sample!
        for jj = jumpinds(:)'
            linkinds = ((jj-1):(jj)) + aa;
            plot(postimevec(linkinds) - postimevec(1),posvec(linkinds - aa + 1),'k-','linewidth',5);
            plot(postimevec(linkinds) - postimevec(1),posvec(linkinds - aa + 1),'-','Color',[40 200 40]/255,'linewidth',2); hold on
        end
    end
    
    % Plot arm boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot([plot_start plot_end],[CPbuff CPbuff],'k--','linewidth',1)
    plot([plot_start plot_end],[outersep outersep],'k--','linewidth',1)
    
    % Plot 1 s grid
    if 0
        for sc = 1:length(plot_start:plot_end)
            plot([plot_start+sc-1 plot_start+sc-1],[0 max(posvec_all)+1],'color',[.5 .5 .5]); hold on
        end
        set(gca,'fontsize',14)
    end
    set(gca,'tickdir','out')
    set(gca,'color','none','fontsize',14)
    

    % Plot LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AX(3) = subplot(7,1,1);
    if plot_DELTA_LFP
        lfp_a = lookup(plot_start,deltatimevec-epstart);
        lfp_b = lookup(plot_end,deltatimevec-epstart);
        plot(deltatimevec(lfp_a:lfp_b)-epstart,LFPSIGN *deltatrace(lfp_a:lfp_b),'-','color',[.75 .75 .9],'linewidth',2); hold on
    end
    if plot_LOWFREQ_LFP
        lfp_a = lookup(plot_start,lowfreqtimevec-epstart);
        lfp_b = lookup(plot_end,lowfreqtimevec-epstart);
        plot(lowfreqtimevec(lfp_a:lfp_b)-epstart,LFPSIGN *lowfreqtrace(lfp_a:lfp_b),'-','color',[.8 .8 .9],'linewidth',2); hold on
    end
    if 1 %plot_THETA_LFP
        lfp_a = lookup(plot_start,thetatimevec-epstart);
        lfp_b = lookup(plot_end,thetatimevec-epstart);
            if any(strcmp(animal_toplot,{'Miles','Conley'}))
                DIVISOR = 0.25;
            else
                DIVISOR = 1;
            end
        plot(thetatimevec(lfp_a:lfp_b)-epstart,DIVISOR*LFPSIGN * thetatrace(lfp_a:lfp_b),'k-','linewidth',1.5); hold on
        ylim([-LFPMAX LFPMAX])
    end
    if plot_LFP
        lfp_a = lookup(plot_start,eegtimevec-epstart);
        lfp_b = lookup(plot_end,eegtimevec-epstart);
        plot(eegtimevec(lfp_a:lfp_b)-epstart,LFPSIGN *eegtrace(lfp_a:lfp_b),'k-','linewidth',2); hold on
    end

    set(gca,'fontsize',14,'color','none')
    set(gca,'xticklabel','')

    % Unit MUA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    AX(4) = subplot(7,1,4);
     muatvec = plot_start:0.001:plot_end;
     muavec = zeros(1,length(muatvec));
     for tet = selected_tets
         filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory
         inds_thresh = any(filedata.params(:,2:5) > P.spikethresh, 2) ;  %  Spikes with min amplitude
         spiketimes = filedata.params(inds_thresh,1)/10000 - epstart;
         onehist = histc(spiketimes,muatvec);
         onehist(end) = 0;
         muavec = muavec + onehist(:)';
     end
     muavec = muavec ;   % in kHz
     msmuakern = 20;  % in ms
     muakern = gaussian(msmuakern,msmuakern*10);
     muavec = smoothvect(muavec,muakern);
     if 0
         figure;
         plot(tvec,muavec,'k-')
     end
    
    area(muatvec,muavec,'facecolor','k');
    set(gca,'xticklabel','')
    set(gca,'tickdir','out');
    set(gca,'fontsize',14,'color','none')
    xlim([plot_start plot_end])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot Speed (linear + angular) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AX(5) = subplot(7,1,5);
    postimevec2 = postimevec - postimevec(1);
    pos_a = lookup(plot_start,postimevec2,-1) - 1;
    pos_b = lookup(plot_end,postimevec2,+1) + 1;
    %  speed
    [postimevec_seg,veltrace_seg] = Areainterp(postimevec2(pos_a:pos_b),veltrace(pos_a:pos_b),[plot_start plot_end]);
    area(postimevec_seg,veltrace_seg,'facecolor',[.72 .72 .72],'linewidth',1);   hold on 
    %  ang speed
    [postimevec_seg,veltrace_seg] = Areainterp(postimevec2(pos_a:pos_b),angvel(pos_a:pos_b),[plot_start plot_end]);
    area(postimevec_seg,ANGVEL_SCALE * veltrace_seg,'facecolor',[.5 .5 .5],'edgecolor','none')
    set(gca,'fontsize',14,'color','none')
       set(gca,'tickdir','out');       
       %set(gca,'children',flipud(get(gca,'children')))
       set(gca,'box','on')
       
       %set(gca,'ticklength',1.5*get(gca,'ticklength'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    linkaxes(AX,'x');

    % plot title
    AX(3) = subplot(7,1,1);
    title(sprintf('%s day %d ep %d',animpref(1:3),d,ep),'fontweight','bold','fontsize',16)

    %keyboard
    %pause
    %close all
    if save_plot
        cd('/opt/data50/kkay/___Fig2')
        %saveas(K,exstr,'pdf')
        export_fig(exstr,'-eps')
    end
    
    end
end
        
    
   
    
   
    cd('/opt/data50/kkay/___Fig2')

% Place map plot: 2D map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_placemap
    F2 = figure('units','normalized','outerposition',[.7 .2 .2 .4]);
    plot_behaviormap(animal_toplot,d,ep,[plot_start plot_end],[])  % plotting function
    title(sprintf('%s day %d ep %d',animalinfo{3}(1:3),d,ep),...
        'fontsize',14,'fontweight','bold')
    set(gca, 'XTick', [],'YTick',[]);
    set(gca,'box','off')
    if save_plot
        PFfilename = [exstr '_map'];
        export_fig(PFfilename,'-r300','-rgb','-tif','-painters')
        %export_fig('Ex_1_placemap','-pdf')
    end
end

break

if save_plot
    close(K)
    close(F2)
end




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


