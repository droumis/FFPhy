% plot locations belonging to one segment on the track

day=7;
epoch=[2];
segment=[1:12];
animals='H2';
dsz = '';
if (day < 10)
    dsz = '0';
end
 dayt = num2str(day);

animaldir='H2_';
datadir = '/data14/jai/';
%posfilename = strcat(datadir,animaldir,'/',animals,'linpos',dsz,dayt,'.mat');
posfilename = strcat(datadir,animaldir,'/',animals,'_',dsz,dayt,'linpos','.mat');
load(posfilename);


segindex=linpos{day}{epoch}.statematrix.segmentIndex;
wellcoord=linpos{day}{epoch}.wellSegmentInfo.wellCoord;
% Load the mat file, needed for barrier information
%sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');
sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'_',dsz,dayt,'.mat');

load(sfilename);

pos=Data{day}{epoch}.Pos.correcteddata;

for segment=segment;

dataind=find(segindex==segment);

datacoord=pos(dataind,2:3);

% distance from point to line

end1=[linpos{day}{epoch}.segmentInfo.segmentCoords(segment,1) linpos{day}{epoch}.segmentInfo.segmentCoords(segment,2)];
end2=[linpos{day}{epoch}.segmentInfo.segmentCoords(segment,3) linpos{day}{epoch}.segmentInfo.segmentCoords(segment,4)];

for i=1:size(datacoord,1);
    point=[datacoord(i,1) datacoord(i,2)];
    d(i) = abs(det([end2-end1;point-end1]))/norm(end2-end1);
end


% figure;
% 
% plot3(datacoord(:,1),datacoord(:,2), d','.r');
% xlim([0 160]);
% ylim([0 120]);
% zlim([0 20]);
% grid on;

figure;
plot(datacoord(:,1),datacoord(:,2),'.r');
hold on;

% get locations of intersections
intersectcoord=unique([linpos{day}{epoch}.segmentInfo.segmentCoords(:,1:2);linpos{day}{epoch}.segmentInfo.segmentCoords(:,3:4)] ,'rows');

wellcoord=intersectcoord;


N=256;
t=(0:N)*2*pi/N;
for i=1:size(wellcoord,1);
    h=wellcoord(i,1); k=wellcoord(i,2); r=5; 
    
    plot( r*cos(t)+h, r*sin(t)+k,'b');
    
    text(wellcoord(i,1),wellcoord(i,2),num2str(i));
    hold on;
    i=i+1;
end


xlim([0 160]);
ylim([0 120]);
title(sprintf('Segmentation for %s for day %s epoch %s segment %s',...
        animals, dayt, num2str(epoch), num2str(segment)));
end
clear all;


    