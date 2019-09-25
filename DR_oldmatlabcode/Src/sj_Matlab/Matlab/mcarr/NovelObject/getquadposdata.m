function [percents quads types]=getquadposdata(day,epoch,minutes)
% function [percents quads types]=getquadposdata(day,epoch,minutes)
% Get the percent of times in each quadrant for a particular epoch and
% day. Also returns a struct where the first column is zeros and ones
% representing for each frame whether the animal was in quadrant one or not
% and the second column the same for quadrant 2, 3rd column q3, 4th column
% q4.  The types return is 4 strings that tell you what the content was of
% columns 1,2,3, and 4 in this particular epoch.
%
% For the minutes input, can either write 'all' to see the whole epoch or
% [a b] to see minutes a through b (including b), i.e. to see just the
% first minute type [1 1], but to see minutes 2, 3, and 4 type [2 4]

%load necessary info
if day==10
    load '/data13/mcarr/NovelObject/Cor/Corpos10.mat';
    load '/data13/mcarr/NovelObject/Cor/Cortask10.mat';
elseif day==11
    load /data13/mcarr/NovelObject/Cor/Corpos11.mat;
    load /data13/mcarr/NovelObject/Cor/Cortask11.mat;
else
    error('invalid day argument');
end

%get x and y values for each quadrant during that epoch
quad1x=task{day}{epoch}.quadrants(1).x;
quad1y=task{day}{epoch}.quadrants(1).y;
quad2x=task{day}{epoch}.quadrants(2).x;
quad2y=task{day}{epoch}.quadrants(2).y;
quad3x=task{day}{epoch}.quadrants(3).x;
quad3y=task{day}{epoch}.quadrants(3).y;
quad4x=task{day}{epoch}.quadrants(4).x;
quad4y=task{day}{epoch}.quadrants(4).y;

%get x and y positions of animal at the times you want
if ischar(minutes)==1    %to see whole epoch    
    xpos=pos{day}{epoch}.data(:,2);
    ypos=pos{day}{epoch}.data(:,3);
else           %to see certain minutes
    framestart=1800*(minutes(1)-1)+1;
    frameend=1800*minutes(2);
    frames=framestart:frameend;
    xpos=pos{day}{epoch}.data(frames,2);
    ypos=pos{day}{epoch}.data(frames,3);
end

%get counts of when in each quadrant
quad1=xpos>=quad1x(1) & xpos<quad1x(2) & ypos>=quad1y(1) & ypos<=quad1y(2);
quad2=xpos>=quad2x(1) & xpos<=quad2x(2) & ypos>=quad2y(1) & ypos<=quad2y(2);
quad3=xpos>=quad3x(1) & xpos<quad3x(2) & ypos>=quad3y(1) & ypos<quad3y(2);
quad4=xpos>=quad4x(1) & xpos<=quad4x(2) & ypos>=quad4y(1) & ypos<quad4y(2);

%check that every time is assigned to one and only one quadrant
total=quad1+quad2+quad3+quad4;
%check that assigned to one
zeros=total==0;
if sum(zeros)~=0
    error('there are times unassigned to quadrants')
end
%check that only one
twos=total>1;
if sum(twos)~=0
    error('there are times assigned to more than one quadrant')
end

%get percent of time in each quadrant (sum of times that the values match
%each quadrant divided by total num of samples)
percentq1=sum(quad1)/size(xpos,1);
percentq2=sum(quad2)/size(xpos,1);
percentq3=sum(quad3)/size(xpos,1);
percentq4=sum(quad4)/size(xpos,1);

%determine each quadrant's type/content
typeq1=task{day}{epoch}.quadrants(1).type;
typeq2=task{day}{epoch}.quadrants(2).type;
typeq3=task{day}{epoch}.quadrants(3).type;
typeq4=task{day}{epoch}.quadrants(4).type;

quads=[quad1 quad2 quad3 quad4];
percents=[percentq1 percentq2 percentq3 percentq4];
types={typeq1 typeq2 typeq3 typeq4};

figure;
h=bar(percents);
hold on;
if ischar(minutes)==1
    tempstring=['Time Spent in Each Quadrant, Day ' num2str(day) ' Epoch ' num2str(epoch) ' All 10 Minutes'];
elseif minutes(1)==minutes(2)
    tempstring=['Time Spent in Each Quadrant, Day ' num2str(day) ' Epoch ' num2str(epoch) ' Minute ' num2str(minutes(1))];
else
    tempstring=['Time Spent in Each Quadrant, Day ' num2str(day) ' Epoch ' num2str(epoch) ' Minutes ' num2str(minutes(1)) ' Through ' num2str(minutes(2))];
end
title(tempstring);
ylabel('Percent of Total Time');
xlabel('Quadrants');
set(gca,'XTickLabel',types);
set(gca,'YLim',[0 1]);
hold off;


