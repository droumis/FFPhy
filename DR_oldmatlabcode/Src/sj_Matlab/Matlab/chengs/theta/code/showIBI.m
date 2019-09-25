function showIBI(plotFigs, dtburst, num)
% function showCoherentSpikes(maxdt, num)
% maxdt= "coherence time length"

global fmaux select dtburst
%plotFigs=  1;
minspikes= 3; % minimum number of spikes or bursts to do meaningful analysis
%dtburst=0.05; % within this time [ms] a spike is considered part of a burst
maxdt=      .6;

format short g
if nargin < 3
    loadFile(fmaux.select);
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c, minspikes,plotFigs,maxdt);
%        disp('any key'); pause
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2),num(3),num(4), minspikes,plotFigs,maxdt);
end


function auxrun(d,e,t,c, minspikes,plotFigs,maxdt)

global fmaux behavdata spikedata select dtburst

loadVar(fmaux.data2dir, 'behavdata', d);
loadVar(fmaux.data2dir, 'spikedata', d);
sd= spikedata{d}{e}{t}{c};
bd= behavdata{d}{e};
if isempty(sd); return; end % no spikes

n= fmaux.currentCell;
%if (n > length(select.x)) | isempty(select.x{n}) ...
%    | (n > length(select.a)) | isempty(select.a{n})
%    return;
%end

xmin=min(bd.linpos);
xmax=max(bd.linpos);

A= select.a{n};
for ia=1:length(A)
    traj= A{ia}.traj;
%for traj=0:3
%    pause
    % select spikes on right traj
    tind= find(bd.traj(sd.index)== traj);

    x= sd.linpos(tind);
    time= sd.time(tind);
    phase= bd.phase(sd.index(tind));
    [C delta]= auxLinCorr(x, phase);
    phase=  mod(phase+delta, 2*pi);
    bd.phase= mod(bd.phase+delta, 2*pi);
    if plotFigs
        if plotFigs==2
            fig= figure(1); clf; set(fig,'Position',[50 400, 1000, 300]);
        else
            fig= figure(1); clf; set(fig,'Position',[700 40, 600, 300]);
        end


        tstring= sprintf('[%d %d %d %d], cellnum= %d, traj= %d, C= %.2f, delta= %.1f\n', d,e,t,c,n,traj, C, delta);
%        plot(sd.linpos(tind), mod(bd.phase(sd.index(tind))+delta, 2*pi), '.'); 
        plot(x, phase, '.'); 
        title(tstring);
        hold on
%        plot(A{ia}.linpos, [0; 0], 'sr');
        hold off
        axis([xmin xmax 0 2*pi]);
        if plotFigs==2
            fig= figure(2); clf; set(fig,'Position',[50 720, 1000, 300]);
        else
            fig= figure(2); clf; set(fig,'Position',[250 300, 1000, 720]);
            fig= figure(3); clf; set(fig,'Position',[250 300, 1000, 720]);
        end
%        keyboard
    end
    
    % find coherent spikes
    if traj==1 | traj==3
        x= -x;
    end
    
    ib= 1;
    nbout= 1;
    xcoh= cell(1,1); pcoh= cell(1,1); tcoh= cell(1,1); icoh= cell(1,1);
    while ib <= length(x)
        xprev= min(x)-1;
        tprev= time(ib);
        xtmp= []; ptmp= []; ttmp= []; itmp= [];
        nburst= 0;
%        while ib <= length(x) & x(ib) > xprev & time(ib)-tprev < maxdt
        while ib <= length(x) & time(ib)-tprev < maxdt
            if nburst== 0 | time(ib)-tprev > dtburst
                xtmp(end+1)= x(ib);
                ptmp(end+1)= phase(ib);
                ttmp(end+1)= time(ib);
                itmp(end+1)= sd.index(tind(ib));
                nburst= nburst+1;
            end
            xprev= x(ib);
            tprev= time(ib);
            ib= ib+1;
        end
        if traj==1 | traj==3 xtmp= -xtmp; end
        if length(xtmp)>= minspikes & ttmp(end)-ttmp(1) > dtburst
            nbout= nbout+1;
            xcoh{nbout}= xtmp;
            pcoh{nbout}= ptmp;
            tcoh{nbout}= ttmp;
            icoh{nbout}= itmp;
        else
            xcoh{1}=[xcoh{1}, xtmp];
            pcoh{1}=[pcoh{1}, ptmp];
            tcoh{1}=[tcoh{1}, ttmp];
            icoh{1}=[icoh{1}, itmp];
        end
    end

    cmap= jet(nbout);
    cmap(1,:)=[1 1 1];
    n= 0;
    waitKey= 1;
    if plotFigs==1
        if(1)
            figure(3)
            for ib=2:nbout
%                if traj==1 | traj==3; tcoh{ib}=max(tcoh{ib})-tcoh{ib}; end
                subplot(ceil(nbout/5),5,ib-1);
                plot(tcoh{ib}-min(tcoh{ib}), pcoh{ib}, 'o', 'MarkerFaceColor',cmap(ib,:));
    %            axis([min(tcoh{ib}) max(tcoh{ib}) 0 2*pi]);
%                axis([0 2 0 2*pi]);
            axis([0 max(tcoh{ib})-min(tcoh{ib})  0 2*pi]);
                hold on
                ind= icoh{ib}(1):icoh{ib}(end);
                plot(bd.time(ind)-min(tcoh{ib}), bd.phase(ind),'.','MarkerSize',1);
                if length(tcoh{1});
                    subplot(ceil(nbout/5),5,nbout);
                    plot(tcoh{1}, pcoh{1}, 'o', 'MarkerFaceColor',cmap(1,:));
%                    axis([0 2 0 2*pi]);
                    axis tight
                end
            end
        end
        figure(2)
        for ib=2:nbout
            subplot(ceil(nbout/5),5,ib-1);
%            plot(xcoh{ib}, pcoh{ib}, '.');
%            [r,deltatmp]= auxLinCorr(xcoh{ib},pcoh{ib});
%            plot(xcoh{ib}, mod(pcoh{ib}+deltatmp,2*pi), 'o', 'MarkerFaceColor',cmap(ib,:));
            plot(xcoh{ib}, pcoh{ib}, 'o', 'MarkerFaceColor',cmap(ib,:));
%            plot(xcoh{ib}, mod(tcoh{ib},0.125), 'o', 'MarkerFaceColor',cmap(ib,:));
            axis([min(xcoh{ib})-1e-3 max(xcoh{ib})+1e-3 0 2*pi]);
%            axis([A{ia}.linpos' 0 2*pi]);
%            axis([xmin xmax 0 2*pi]);
%            axis([min(xcoh{ib}) max(xcoh{ib}) 0 0.125]);
            hold on
            ind= icoh{ib}(1):icoh{ib}(end);
%            plot(bd.linpos(ind), mod(bd.phase(ind)+deltatmp,2*pi),'.','MarkerSize',1);
            plot(bd.linpos(ind), bd.phase(ind),'.','MarkerSize',1);
%            keyboard
            if length(xcoh{1});
                subplot(ceil(nbout/5),5,nbout);
%                 plot(xcoh{1}, mod(pcoh{1}+deltatmp,2*pi), 'o', 'MarkerFaceColor',cmap(1,:));
            plot(xcoh{1}, pcoh{1}, 'o', 'MarkerFaceColor',cmap(1,:));
            axis([xmin xmax 0 2*pi]);
            end
        end
        orient landscape
%        print -dps -P811
    elseif plotFigs==2
        for ib=2:nbout
            hold on
%            plot(xcoh{ib}, pcoh{ib}, 'o', 'MarkerFaceColor',cmap(ib,:));
%            axis([xmin xmax 0 2*pi]);
            plot(tcoh{ib}-min(tcoh{ib}), mod(pcoh{ib}-pcoh{ib}(1)-1e-7,2*pi), 'o', 'MarkerFaceColor',cmap(ib,:));
            axis([0 2 0 2*pi]);
            if waitKey
                r=input('(s)kip to next traj: ','s');
                if r=='s'; waitKey= 0; end
            end
        end
        if length(xcoh{1});
%            plot(xcoh{1}, pcoh{1}, 'o', 'MarkerFaceColor',cmap(1,:));
            plot(tcoh{1}, pcoh{1}, 'o', 'MarkerFaceColor',cmap(1,:));
            axis tight
%            axis([0 2 0 2*pi]);
        end
    else
        for ib=1:nbout
%            [tcoh{ib}'-min(tcoh{ib}) xcoh{ib}' pcoh{ib}']
            fprintf('%2d \t%g \t%g \t%g\n', ib, tcoh{ib}(1), xcoh{ib}(1),...
                pcoh{ib}(1));
        end
    end
    disp('done'); pause
end

