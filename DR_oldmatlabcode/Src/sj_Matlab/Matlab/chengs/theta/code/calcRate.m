function data=calcRate(num)
%function data= calcRate(num)
%  
%  Calculated mean firing rate as a function of linear distance.
%    data   column vector with center of spatial bin and rate for each of the
%           four trajectories.

global fmaux
if nargin < 1
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        rate{d}{e}{t}{c}.fields= 'linpos rate1 rate2 rate3 rate4';
        rate{d}{e}{t}{c}.data= auxrun(d,e,t,c);
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2),num(3),num(4));
end

save([fmaux.data2dir '/rate'], 'rate');


function data= auxrun(d,e,t,c)

xbin= .5; % [cm]
occmin= 0.050; % [sec]
gstd= 4;
g= gaussian(gstd, 6*gstd+1);
pad= 1; % [cm]

global fmaux spikedata behavdata

loadVar(fmaux.data2dir, 'spikedata', d);
loadVar(fmaux.data2dir, 'behavdata', d);

bd= behavdata{d}{e};
sd= spikedata{d}{e}{t}{c};

xmin= min(bd.linpos)-pad;
xmax= max(bd.linpos)+pad;
nbin= round((xmax-xmin)/xbin);
x= linspace(xmin, xmax, nbin+1)';

dt= mean(diff(bd.time));
data= zeros(nbin, 5);
data(:,1)= x(1:end-1)+(x(2)-x(1))/2;

for traj=0:3
    data(:,traj+2)= zeros(nbin,1);
    ind= find(bd.traj== traj); 
    if isempty(ind); continue; end
    occ= histc(bd.linpos(ind), x);
    if(sum(occ)~=length(ind)) error('linpos values outside [xmin xmax]'); end
    valid= bd.traj(sd.index)== traj;
    if sum(valid)<2; continue; end
    nsp= histc(sd.linpos(valid),x); 
    if(sum(nsp)~=sum(valid)) error('spike linpos values outside [xmin xmax]'); end
    occ= occ*dt;
    rate= zeros(nbin,1);
    nonzero= find(occ(1:end-1)>occmin);
    if isempty(nonzero); continue; end
    rate(nonzero)= nsp(nonzero)./occ(nonzero);
    data(:,traj+2)=filtfilt(g,1,rate);
%    plot(x(1:end-1), [rate data(:,traj+2)]);
%    pause 
end

if 0 
    figure(1)
%        subplot(2,1,2);
%        plot(sd.linpos(valid), bd.phase(sd.index(valid)), '.'); 
%        axis([0 160, 0 2*pi]);
%        subplot(2,1,1);
%    tstring= sprintf('[%d %d %d %d], traj= %d\n',d,e,t,c,traj);
%    plot(data(:,1), [rate, data(:,traj+2)]);
%    axis([0 160, 0 max(max([rate, data(:,traj+2)]))]);

    tstring= sprintf('[%d %d %d %d]\n',d,e,t,c);
    plot(data(:,1), data(:,2:5));
    axis([0 160, 0 max(max(data(:,2:5)))]);
    legend({'rate', 'vel'});
    title(tstring);
    pause
end
