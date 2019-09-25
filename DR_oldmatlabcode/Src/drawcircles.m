
animal='N2';
day= 1;
epoch= 4;


dsz = '';
if (day < 10)
    dsz = '0';
end
dayt = num2str(day);
datadir = '/data14/jai/';
sfilename = strcat(datadir,animal,'_/',animal,'data',dsz,dayt,'.mat');
posfilename = strcat(datadir,animal,'_/',animal,'task',dsz,dayt,'.mat');
load(sfilename);
load(posfilename);

% wellcoord=linpos{1,day}{1,epoch}.wellSegmentInfo.wellCoord;
%wellcoord=unique([linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,1:2);linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,3:4)],'rows');
%wellcoord=linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,1:2);
task{1,day}{1,epoch}.linearcoord;

allintersect=[];

for i = 1:length(task{1,day}{1,epoch}.linearcoord)
    coordinates = task{1,day}{1,epoch}.linearcoord{i}(:,:,1);
    allintersect=[allintersect;coordinates];
end

wellcoord=unique(allintersect,'rows');




N=256;
fid=figure;
set(fid,'OuterPosition',[1 1 1280 1024])
for e=epoch;

%plot(data{1,day}{1,e}.Pos.rawpos(:,2)*data{1,day}{1,e}.Pos.cmperpixel/2,data{1,day}{1,e}.Pos.rawpos(:,3)*data{1,day}{1,e}.Pos.cmperpixel/2,'.','Color',[0.2 0.2 0.2]);

plot(data{1,day}{1,e}.Pos.rawpos(:,2),data{1,day}{1,e}.Pos.rawpos(:,3),'.','Color',[0.2 0.2 0.2]);

hold on;
end

% for e=epoch
% 
% plot(data{1,day}{1,e}.Pos.correcteddata(:,2),data{1,day}{1,e}.Pos.correcteddata(:,3),'.','Color',[0.5 0.5 0.5]);
% hold on;
% end


%plot(data{1,day}{1,epoch}.Pos.correcteddata(:,2)*data{1,day}{1,e}.Pos.cmperpixel,data{1,day}{1,epoch}.Pos.correcteddata(:,3)*data{1,day}{1,e}.Pos.cmperpixel,'.r');
plot(data{1,day}{1,epoch}.Pos.correcteddata(:,2),data{1,day}{1,epoch}.Pos.correcteddata(:,3),'.r');

hold on;

t=(0:N)*2*pi/N;
for i=1:size(wellcoord,1);
    h=wellcoord(i,1); k=wellcoord(i,2); r=15; 
    
    plot( r*cos(t)+h, r*sin(t)+k,'b');
    
    %text(wellcoord(i,1),wellcoord(i,2),num2str(i));
    hold on;
    i=i+1;
end
xlim([0 320]);
ylim([0 240]);
title(sprintf('Day %d epoch %d',day,epoch));

    