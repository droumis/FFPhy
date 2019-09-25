function [riprates types]=getripplerate(day,epoch,minutes)
% function out=getripplerate(day,epoch,minutes)
% Get ripples/second in a particular epoch in each quadrant.
%
% For the minutes input, can either write 'all' to see the whole epoch or
% [a b] to see minutes a through b (including b), i.e. to see just the
% first minute type [1 1], but to see minutes 2, 3, and 4 type [2 4]


% load necessary info
load '/data13/mcarr/NovelObject/Cor/Corcellinfo.mat';
if day==10
    load '/data13/mcarr/NovelObject/Cor/Corpos10.mat';
    load '/data13/mcarr/NovelObject/Cor/Cortask10.mat';
    load '/data13/mcarr/NovelObject/Cor/Corripples10.mat';
elseif day==11
    load /data13/mcarr/NovelObject/Cor/Corpos11.mat;
    load /data13/mcarr/NovelObject/Cor/Cortask11.mat;
    load /data13/mcarr/NovelObject/Cor/Corripples11.mat;
else
    error('invalid day argument');
end

% get the times at which ripples occurred
rippletimes = getripples([day epoch],ripples,cellinfo,'cellfilter', 'isequal($area,''CA1'')');
% get what index each ripple is in the position data
rippleposinds = lookup(rippletimes(:,1), pos{day}{epoch}.data(:,1));

% if dont want whole epoch, modify rippleposinds to only include ripples
% that occurred during the requested minutes.
if ischar(minutes)==0  % when only want certain minutes
    framestart=1800*(minutes(1)-1)+1;
    frameend=1800*minutes(2);
    j=1;
    for i=1:size(rippleposinds)
        if rippleposinds(i)>=framestart && rippleposinds(i)<=frameend
            temp(j,1)=rippleposinds(i);
            j=j+1;
        end
    end
    rippleposinds=temp;
end

% get list of when rat was in each quadrant, need to do all minutes here
% because otherwise the rippleposinds wont line up with the correct place
% in quads when starting from minute > 1.
[~, quads t] = getquadposdata(day,epoch,'all'); 
% close the figure because dont want to see position for whole epoch
close;

% get list of ripples that happen in each quadrant
q1ripples=quads(rippleposinds,1)==1;
q2ripples=quads(rippleposinds,2)==1;
q3ripples=quads(rippleposinds,3)==1;
q4ripples=quads(rippleposinds,4)==1;

%check that every ripple is assigned to one and only one quadrant
total=q1ripples+q2ripples+q3ripples+q4ripples;
%check that assigned to one
nones=total==0;
if sum(nones)~=0
    error('there are ripples unassigned to quadrants')
end
%check that only one
twos=total>1;
if sum(twos)~=0
    error('there are ripples assigned to more than one quadrant')
end

% get ripple rate-- total number of ripples in a quadrant divided by the
% total number of frames the animal was in that quadrant, times 30
% frames/second.
q1riprate=30*sum(q1ripples)/sum(quads(:,1));
q2riprate=30*sum(q2ripples)/sum(quads(:,2));
q3riprate=30*sum(q3ripples)/sum(quads(:,3));
q4riprate=30*sum(q4ripples)/sum(quads(:,4));

riprates=[q1riprate q2riprate q3riprate q4riprate];

%determine each quadrant's type/content
typeq1=task{day}{epoch}.quadrants(1).type;
typeq2=task{day}{epoch}.quadrants(2).type;
typeq3=task{day}{epoch}.quadrants(3).type;
typeq4=task{day}{epoch}.quadrants(4).type;

types={typeq1 typeq2 typeq3 typeq4};

% do this to get graph of position in just the requested minutes
getquadposdata(day,epoch,minutes);

figure;
h=bar(riprates);
hold on;
if ischar(minutes)==1
    tempstring=['Ripple Rate in Each Quadrant, Day ' num2str(day) ' Epoch ' num2str(epoch) ' All 10 Minutes'];
elseif minutes(1)==minutes(2)
    tempstring=['Ripple Rate in Each Quadrant, Day ' num2str(day) ' Epoch ' num2str(epoch) ' Minute ' num2str(minutes(1))];
else
    tempstring=['Ripple Rate in Each Quadrant, Day ' num2str(day) ' Epoch ' num2str(epoch) ' Minutes ' num2str(minutes(1)) ' Through ' num2str(minutes(2))];
end
title(tempstring);
ylabel('Ripples/Second');
xlabel('Quadrants');
set(gca,'XTickLabel',types);
hold off;
