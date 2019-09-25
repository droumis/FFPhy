function showBehav(data,fig,Tmax)

if nargin < 2
    fig= 1;
end
if nargin < 3
    Tmax= length(data.timearray);
else
    if Tmax > length(data.timearray)
	warning('requested more data than available in file');
    end
end

if Tmax ~= length(data.timearray)
    data.timearray= data.timearray(1:Tmax);
    data.posarray= data.posarray(1:Tmax);
    data.phasearray= data.phasearray(1:Tmax);
    data.fieldID= data.fieldID(1:Tmax);
end

x=linspace(min(data.posarray),max(data.posarray),15);
y=linspace(min(data.phasearray),max(data.phasearray),15);

dx=x(2)-x(1);
dy=y(2)-y(1);
dT= (data.timearray(2)-data.timearray(1));
occup= hist2d([data.posarray, data.phasearray], x, y)*dT;

figure(fig); clf;
subplot(3,2,1);
plot(data.timearray, data.fieldID);
title('trajectory');

subplot(3,2,3);
plot(data.timearray, data.posarray);
title('position');

subplot(3,2,5);
plot(data.timearray, data.phasearray);
title('theta phase');

subplot(3,2,2);
%figure(fig+1); clf;
imagesc(x,y,(occup)');
axis xy
colorbar; title('occupancy');
