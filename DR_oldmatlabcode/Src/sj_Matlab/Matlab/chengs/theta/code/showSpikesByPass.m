function showSpikesByPass(plotFigs, num)

minspikes= 4; %minimum number of spikes to do meaningful analysis
%plotFigs=  2;


global fmaux select 
if nargin < 2
    loadFile(fmaux.select);
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c, minspikes,plotFigs);
        disp('press key for next cell'); pause
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2),num(3),num(4), minspikes,plotFigs);
end


function auxrun(d,e,t,c, minspikes,plotFigs)

global fmaux behavdata spikedata select 

loadVar(fmaux.data2dir, 'behavdata', d);
loadVar(fmaux.data2dir, 'spikedata', d);
sd= spikedata{d}{e}{t}{c};
bd= behavdata{d}{e};
if isempty(sd); return; end % no spikes

trajs= bd.traj(sd.index);
n= fmaux.currentCell;
if (n > length(select.x)) | isempty(select.x{n}) ...
    | (n > length(select.a)) | isempty(select.a{n})
    return;
end
xmin= min(bd.linpos);
xmax= max(bd.linpos);

A= select.a{n};
for ia=1:length(A)
%    pause
    % select timesteps on right traj
    itraj= A{ia}.traj;
    tind= find(trajs== itraj);

    % narrow down selection to timesteps when rat was within placefield
%    x= sd.linpos(tind);
%    xind= find(x > A{ia}.linpos(1) & x < A{ia}.linpos(2));
%    xtind= tind(xind);
    xtind= tind; %@@ don't narrow selection down to placefield
    x= sd.linpos(xtind);
    phase= bd.phase(sd.index(xtind));
    [C delta]= auxLinCorr(x, phase);
    phase= mod(phase+delta, 2*pi);

    if plotFigs
        figure(1); clf
        tstring= sprintf('[%d %d %d %d], cellnum= %d, traj= %d, C= %.2f, delta= %.1f\n', d,e,t,c,n,itraj, C, delta);
%        plot(sd.linpos(tind), mod(bd.phase(sd.index(tind))+delta, 2*pi), '.'); 
        plot(x, phase, '.'); 
        title(tstring);
        hold on
%        plot(A{ia}.linpos, [0; 0], 'sr');
        hold off
        axis tight;
        figure(2); clf
    end

    % find analysis intervals
    endtimes= select.x{n}{ia}.time;
    nint= length(endtimes);
    starttimes= [bd.time(1) endtimes(1:nint-1)];

    cmap= jet(nint);

    waitkey=1;
    for iint=1:nint
        % select the spiketimes within one time period (e.g. pass)
        ttind= find(sd.time(xtind) < endtimes(iint) & sd.time(xtind) > starttimes(iint));
        xttind= xtind(ttind);
        x= sd.linpos(xttind);
        phase= bd.phase(sd.index(xttind));
        time= bd.time(sd.index(xttind));

        nspikes= length(x);
        if  nspikes >= minspikes 
%                [C deltatmp]= auxLinCorr(x, phase);
%                phase= mod(phase+deltatmp, 2*pi);
                if plotFigs==1
                    subplot(ceil(nint/4),4,iint);
                    hold off
                    plot(x, phase, 'o', 'MarkerSize', 4, 'MarkerEdgeColor','k','MarkerFaceColor', cmap(iint,:)); 
                    axis([min(x) max(x) 0 2*pi ]);
                elseif plotFigs==2
                    hold on
                    plot(x, phase, 'o', 'MarkerSize', 4, 'MarkerEdgeColor','k','MarkerFaceColor', cmap(iint,:)); 
%                tstring= sprintf('[%d %d %d %d], cellnum= %d, traj= %d, C=
%                %.2f, delta= %.1f\n', d,e,t,c,n,itraj, C, deltatmp);
%                title(tstring);
%                    plot(A{ia}.linpos, [0; 0], 'sr');
                    axis([xmin xmax 0 2*pi ]);
                    if waitkey
                        r= input('(s)kip to next traj','s');
                        if r=='s'; waitkey=0; end
                    end
                else  % text output
                    [time-min(time),x, phase]
                    xtmp= x(1); ptmp= phase(1); ttmp= time(1);
                    for z=2:nspikes
                        if time(z)-time(z-1) > 0.02
                            xtmp(end+1)= x(z);
                            ptmp(end+1)= phase(z);
                            ttmp(end+1)= time(z);
                        end
                    end
                    [ttmp'-min(ttmp), xtmp' ptmp' [0;diff(ttmp')] [0;diff(ptmp')]]
                    [corrcoef(xtmp,ptmp), corrcoef(ttmp,ptmp)]
                    [ptmp'\[xtmp' ones(length(xtmp),1)], ...
                    ptmp'\[ttmp' ones(length(xtmp),1)] ]
                    r= input('k for keyboard: ','s');
                    if r=='k'; keyboard; end
                    if r=='n'; break; end
                    if r=='s'; return; end
                end
        end
    end
    disp('next traj'); pause
end

% shift phases such that |C| becomes maximal
function [Cmax, phimax]= auxLinCorr(x,phase)
stepsize= 0.1;

validind= find(isfinite(phase) & isfinite(x));
x= x(validind);
phase= phase(validind);
if(length(unique(x))<=1 | length(unique(phase))<=1); 
    Cmax= nan;
    phimax= nan;
    return;
end
Cmax= corrcoef(x,phase);
Cmax= Cmax(1,2);
phimax= 0;
for phi=stepsize:stepsize:2*pi
    C= corrcoef(x,mod(phase+phi,2*pi));
    C= C(1,2);
    if abs(C) > abs(Cmax) | (isnan(Cmax) & ~isnan(C))
        Cmax= C;
        phimax= phi;
    end
end
