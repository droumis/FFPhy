function [psth,x,y]= showHist(data,xbin,ybin,fig,Tmax)

% [psth,x,y]= showHist(data,xbin,ybin,fig,Tmax)
%

if nargin < 2
    xbin= 10;
end
if nargin < 3
    ybin= 10;
end
if nargin < 4
    fig= 1;
end
if nargin < 5
    Tmax= length(data.timearray);
else
    if Tmax > length(data.timearray)
	warning('requested more data than available in file');
    end
    data= selectData(data,Tmax);
end

pos= data.posarray(data.ispikes);
phase= data.phasearray(data.ispikes);

x=linspace(min(data.posarray),max(data.posarray),xbin);
y=linspace(min(data.phasearray),max(data.phasearray),ybin);
m2d=hist2d([pos,phase],x,y);

dx=x(2)-x(1);
dy=y(2)-y(1);
% time spend in one cell of histogram
dT= (data.timearray(2)-data.timearray(1));
%dT= dx/(data.posarray(2)-data.posarray(1))*(data.timearray(2)-data.timearray(1))*10;
occup= hist2d([data.posarray, data.phasearray], x, y)*dT;
psth= m2d./occup;

figure(fig); clf;
imagesc(x,y,(psth)');
%imagesc(x,y,m2d'/(dT*dy));
colorbar; title('PSTH');
ylabel('\theta phase [bin]'); xlabel('position [bin]');
axis xy


% figure(fig+1); clf;
% subplot(3,1,1);
% plot(data.timearray, data.posarray);
% 
% subplot(3,1,2);
% plot(data.timearray, data.phasearray);
% 
% subplot(3,1,3);
% hist(data.spiketimes,25);
% [N,X]= hist(data.spiketimes,25);
% N=N/(X(2)-X(1));
% bar(X,N);
%disp('mean firing rate= '); disp(mean(N));

% figure(fig+2); clf;
% imagesc(x,y,(m2d)');
% axis xy
% colorbar; title('raw spike count');
% 
% figure(fig+3); clf;
% imagesc(x,y,(occup)');
% axis xy
% colorbar; title('occupancy');
