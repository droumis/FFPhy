function showIBI(plotFigs, dtburst, num)
% function showCoherentSpikes(maxdt, num)
% maxdt= "coherence time length"

global fmaux select ibi var1 var2
%plotFigs=  1;
minspikes= 3; % minimum number of spikes or bursts to do meaningful analysis
%dtburst=0.05; % within this time [ms] a spike is considered part of a burst
maxdt=      1;

ibi= [];
var1= [];
var2= [];
format short g
if nargin < 3
    loadFile(fmaux.select);
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c, minspikes,plotFigs,maxdt,dtburst);
%        disp('any key'); pause
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2),num(3),num(4), minspikes,plotFigs,maxdt,dtburst);
end

%auxShow([fmaux.prefix ' ibi'], ibi, sqrt(2), plotFigs);
auxShow([fmaux.prefix ' var1'], var1, sqrt(2), plotFigs);
auxShow([fmaux.prefix ' var2'], var2, 1, plotFigs);
auxShowXC([fmaux.prefix ' var2 vs var1'], var1, var2, plotFigs);
%keyboard

function auxShowXC(name, d1, d2, plotFigs)
if isempty(d1)|isempty(d2); disp('data empty'); return; end
C=corrcoef(d1,d2);
fprintf(1,'%s: lin correl coeff= %.3f\n', name, C(1,2));
switch plotFigs
case 0,
case 1,
    figure
    plot(d1,d2,'.');
    tstring= sprintf('%s',name);
    title(tstring)
    axis tight
otherwise
    warning(['plotFig has unknown value ' num2str(plotFig)]);
end

function auxShow(name, data, factor, plotFigs)

if isempty(data); disp('data empty'); return; end
sd= std(data);
md= mean(data);
fprintf(1,'%s: n= %d, mean= %.2f, std= %.2f(-> %.2f)\n',name,length(data),md,sd, sd/factor);
b=[-inf, linspace(-2*sd,2*sd,20), inf];
h= histc(data-md,b);
switch plotFigs
case 0,
    h
case 1,
    figure
    binwidth= b(3)-b(2);
    bleft=b(2)-2*binwidth;
    bright=b(end-1)+(b(3)-b(2));
    bar([bleft b(2:end-2) bright]-binwidth/2,h(1:end-1));
    tstring= sprintf('%s: mean= %.2f, std= %.2f(-> %.2f)\n',name,md, sd, sd/factor);
    title(tstring)
    axis tight
%    hh = findobj(gca,'Type','patch');
%    set(hh,'FaceColor','k','EdgeColor','w')
otherwise
    warning(['plotFig has unknown value ' num2str(plotFig)]);
end

function auxrun(d,e,t,c, minspikes,plotFigs,maxdt,dtburst)

global fmaux behavdata spikedata select var1 ibi var2

loadVar(fmaux.data2dir, 'behavdata', d);
loadVar(fmaux.data2dir, 'spikedata', d);
sd= spikedata{d}{e}{t}{c};
bd= behavdata{d}{e};
if isempty(sd); return; end % no spikes

n= fmaux.currentCell;
A= select.a{n};
for ia=1:length(A)
    traj= A{ia}.traj;
%for traj=0:3
%    pause
    % select spikes on right traj
    tind= find(bd.traj(sd.index)== traj);
    if length(tind) < minspikes; continue; end

    x= sd.linpos(tind);
    time= sd.time(tind);
    phase= bd.phase(sd.index(tind));
    [C delta]= auxLinCorr(x, mod(phase,2*pi));
    phase= phase+delta;

    ib= 1;
    xprev= min(x)-1;
    tprev= time(ib);
    xtmp= []; ptmp= []; ttmp= []; itmp= [];
    nburst= 0;
    while ib <= length(x)
%        if nburst== 0 | time(ib)-tprev > dtburst & x(ib) > A{ia}.linpos(1) & x(ib) < A{ia}.linpos(2)
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

    d= diff(ptmp);
    ncycles= [diff(floor(ptmp/(2*pi))), 1e7];
    dx= diff(xtmp);
    dt= 1000*diff(ttmp);
    valind= find(ncycles <= 1);
    vald= dt(valind);
    ibi= [ibi, vald];

% count number of bursts in sequence
%%%%%%%%%%%%%%%%%%%%
    nburstseq= zeros(1,length(ttmp));
%    d21clean= [];
    invalid= [0, find(ncycles >= 2), length(d)+1];
    for i= 1:length(invalid)-1
        range= invalid(i)+1:invalid(i+1)-1;
        if ~isempty(range) 
            nburstseq(range)= length(range)+1;
%               d21clean= [d21clean diff(dt(range))];
       end
    end

    valind3= find(nburstseq >= 3);
%    valind2= find(d(1:end-1) > 0*pi & d(1:end-1) < 3*pi & d(2:end) > 0*pi & d(2:end) < 3*pi);

% select only first burst in cycle
%%%%%%%%%%%%%%%%%%%%
    first=[]; i= 1;
    while i<=length(ttmp)
        first(end+1)= i;
        while i < length(ttmp) & ncycles(i+1)==0; i= i+1; end
        i= i+1;
    end
    ncf= diff(floor(ptmp(first)/(2*pi)));
    dtf= diff(1000*ttmp(first));
    dxf= diff(xtmp(first));
    pf= ptmp(first);

% get slope directly
    izerodx= find(abs(dxf)<0.1);
    dxtmp= dxf;
    dxtmp(izerodx)= nan;
    m= diff(mod(pf,2*pi))./dxtmp;

% velocity, hack: linpos was discretized, velocity might not be const.
    vel= 1000*dxf./dtf;

%    val1= dt(valind);
%    val1= ncycles(1:end-1);
%    val2= dtf(find(ncf<= 1));
%    val2= dt;
    val1= m;
    val2= dtf;

%    keyboard
    igood= find(isfinite(val1) & isfinite(val2));
    val1= val1(igood);
    val2= val2(igood);

%    auxShow('ibi', vald, sqrt(2), plotFigs);
%    auxShow('var1', val1, sqrt(2), plotFigs);
%    auxShow('var2', val2, sqrt(2), plotFigs);
%    auxShowXC([fmaux.prefix ' var2 vs var1'], val1, val2, plotFigs);

    var1=[var1, val1];
    var2=[var2, val2];
%    keyboard
end
