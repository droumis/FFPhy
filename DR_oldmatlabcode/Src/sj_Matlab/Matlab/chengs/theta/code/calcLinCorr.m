function calcLinCorr(num)

global stats fmaux select 

minspikes= 10; %minimum number of spikes to do meaningful analysis

fname= ['stats-' fmaux.selectid '.mat'];
if exist(fname) > 0
    load(fname);
else
    stats= [];
end
loadFile(fmaux.select);

if nargin < 1
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c, minspikes);
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2),num(3),num(4), minspikes);
end

save(fname, 'stats');

function auxrun(d,e,t,c, minspikes)

global fmaux behavdata spikedata stats select task

loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, true);
loadVar(fmaux.data2dir, 'behavdata', d);
loadVar(fmaux.data2dir, 'spikedata', d);
sd= spikedata{d}{e}{t}{c};
bd= behavdata{d}{e};
if isempty(sd); return; end % no spikes

trajs= bd.traj(sd.index);
n= fmaux.currentCell;
if (n > length(select.x)) | isempty(select.x{n}) ...
    | (n > length(select.a)) | isempty(select.a{n})
    stats{d}{e}{t}{c}.r= [];
    stats{d}{e}{t}{c}.nspikes= [];
    stats{d}{e}{t}{c}.delta= [];
    return;
end

A= select.a{n};
for ia=1:length(A)
%    pause
    % select timesteps on right traj
    itraj= A{ia}.traj;
    tind= find(trajs== itraj);

    % narrow down selection to timesteps when rat was within placefield
    x= sd.linpos(tind);
    xind= find(x > A{ia}.linpos(1) & x < A{ia}.linpos(2));
    xtind= tind(xind);
    x= sd.linpos(xtind);
    phase= bd.phase(sd.index(xtind));
    [C delta]= auxLinCorr(x, phase);
%delta=5.4; %@@

    if 0 % for debugging
        tstring= sprintf('[%d %d %d %d], cellnum= %d, traj= %d, C= %.2f, delta= %.1f\n', d,e,t,c,n,itraj, C, delta);
%        plot(sd.linpos(tind), mod(bd.phase(sd.index(tind))+delta, 2*pi), '.'); 
        plot(x, mod(phase+delta, 2*pi), '.'); 
        title(tstring);
        hold on
        plot(A{ia}.linpos, [0; 0], 'sr');
        hold off
        axis tight;
        pause
    end

    % find analysis intervals
    endtimes= select.x{n}{ia}.time;
    nint= length(endtimes);
    starttimes= [bd.time(1) endtimes(1:nint-1)];

    stats{d}{e}{t}{c}.r{ia}= zeros(nint,1);
    stats{d}{e}{t}{c}.nspikes{ia}= zeros(nint,1);
    stats{d}{e}{t}{c}.delta{ia}= zeros(nint,1);

for iint=1:nint
        % select the spiketimes within one time period (e.g. pass)
        ttind= find(sd.time(xtind) < endtimes(iint) & sd.time(xtind) > starttimes(iint));
        xttind= xtind(ttind);
        x= sd.linpos(xttind);
        phase= bd.phase(sd.index(xttind));
        time= bd.time(sd.index(xttind));

        nspikes= length(x);
        C= NaN;
        if  nspikes >= minspikes 
            C= corrcoef(x,mod(phase+delta,2*pi));
            C= C(1,2);
%            [C delta]= auxLinCorr(x, phase);
            if(itraj==1 | itraj==3); C= -C; end

            if 0 % for debugging
                tstring= sprintf('[%d %d %d %d], cellnum= %d, traj= %d, C= %.2f, delta= %.1f\n', d,e,t,c,n,itraj, C, delta);
        %        plot(sd.linpos(tind), mod(bd.phase(sd.index(tind))+delta, 2*pi), '.'); 
                plot(x, mod(phase+delta, 2*pi), '.'); 
                title(tstring);
                hold on
                plot(A{ia}.linpos, [0; 0], 'sr');
                hold off
                axis tight;
                pause
            end
            phase= mod(phase+delta, 2*pi);
            [time-min(time),x, phase]
            xtmp= x(1); ptmp= phase(1); ttmp= time(1);
            for z=2:nspikes
                if time(z)-time(z-1) > 0.02
                    xtmp(end+1)= x(z);
                    ptmp(end+1)= phase(z);
                    ttmp(end+1)= time(z);
                end
            end
            [xtmp' ptmp']
            corrcoef(xtmp,ptmp)
            ptmp'\[xtmp' ones(length(xtmp),1)]
            keyboard

        end

        stats{d}{e}{t}{c}.LinCorr{ia}(iint)= C;
        stats{d}{e}{t}{c}.nspikes{ia}(iint)= nspikes;
        stats{d}{e}{t}{c}.delta{ia}(iint)= delta;
    end
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
