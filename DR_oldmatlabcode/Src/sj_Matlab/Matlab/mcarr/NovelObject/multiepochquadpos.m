function out=multiepochquadpos(minutes)
% combines the behavioral data (when in each quadrant) across multiple
% epochs, correcting for the rotation of objects.
%
% For the minutes input, can either write 'all' to see the whole epoch or
% [a b] to see minutes a through b (including b), i.e. to see just the
% first minute type [1 1], but to see minutes 2, 3, and 4 type [2 4]

% would ultimately want to be able to feed particular days and epochs, but
% for now, just hard coded since only 4. and would want to be able to feed
% a variable number.
epochs=[10 4;10 8;11 4;11 8];

% extract the desired day and epoch numbers
dayA=epochs(1,1);
epochA=epochs(1,2);
dayB=epochs(2,1);
epochB=epochs(2,2);
dayC=epochs(3,1);
epochC=epochs(3,2);
dayD=epochs(4,1);
epochD=epochs(4,2);
%numberepochs=size(epochs,1);

% get quadrant position data for each desired epoch
[percentsA, ~, typesA]=getquadposdata(dayA,epochA,minutes);
[percentsB, ~, typesB]=getquadposdata(dayB,epochB,minutes);
[percentsC, ~, typesC]=getquadposdata(dayC,epochC,minutes);
[percentsD, ~, typesD]=getquadposdata(dayD,epochD,minutes);
close all; %to close the figures for each epoch, can always comment this out if want to see them

% find which quadrant was the NovObj quad in each epoch
whereNOa=strcmp('NovObj',typesA); %where is the NovObj in epoch A
whereNOb=strcmp('NovObj',typesB);
whereNOc=strcmp('NovObj',typesC);
whereNOd=strcmp('NovObj',typesD);
% find the percent of times in NovObj quad for each epoch and put those
% together. 
NovObjPercents=[0 0 0 0];
NovObjPercents(1)=percentsA(1,whereNOa);
NovObjPercents(2)=percentsB(1,whereNOb);
NovObjPercents(3)=percentsC(1,whereNOc);
NovObjPercents(4)=percentsD(1,whereNOd);
NovObjAvg=nanmean(NovObjPercents);


% find which quadrant was the FamObj quad in each epoch
whereFOa=strcmp('FamObj',typesA); %where is the FamObj in epoch A
whereFOb=strcmp('FamObj',typesB);
whereFOc=strcmp('FamObj',typesC);
whereFOd=strcmp('FamObj',typesD);
% find the percent of times in FamObj quad for each epoch and put those
% together. 
FamObjPercents=[0 0 0 0];
FamObjPercents(1)=percentsA(1,whereFOa);
FamObjPercents(2)=percentsB(1,whereFOb);
FamObjPercents(3)=percentsC(1,whereFOc);
FamObjPercents(4)=percentsD(1,whereFOd);
FamObjAvg=nanmean(FamObjPercents);


% find which quadrant was the NovEmpty quad in each epoch
whereNEa=strcmp('NovEmpty',typesA); %where is the NovEmpty in epoch A
whereNEb=strcmp('NovEmpty',typesB);
whereNEc=strcmp('NovEmpty',typesC);
whereNEd=strcmp('NovEmpty',typesD);
% find the percent of times in NovEmpty quad for each epoch and put those
% together. 
NovEmptyPercents=[0 0 0 0];
NovEmptyPercents(1)=percentsA(1,whereNEa);
NovEmptyPercents(2)=percentsB(1,whereNEb);
NovEmptyPercents(3)=percentsC(1,whereNEc);
NovEmptyPercents(4)=percentsD(1,whereNEd);
NovEmptyAvg=nanmean(NovEmptyPercents);


% find which quadrant was the FamEmpty quad in each epoch
whereFEa=strcmp('FamEmpty',typesA); %where is the FamEmpty in epoch A
whereFEb=strcmp('FamEmpty',typesB);
whereFEc=strcmp('FamEmpty',typesC);
whereFEd=strcmp('FamEmpty',typesD);
% find the percent of times in FamEmpty quad for each epoch and put those
% together. 
FamEmptyPercents=[0 0 0 0];
FamEmptyPercents(1)=percentsA(1,whereFEa);
FamEmptyPercents(2)=percentsB(1,whereFEb);
FamEmptyPercents(3)=percentsC(1,whereFEc);
FamEmptyPercents(4)=percentsD(1,whereFEd);
FamEmptyAvg=nanmean(FamEmptyPercents);


out=[NovObjAvg FamObjAvg NovEmptyAvg FamEmptyAvg];

figure;
h=bar(out);
hold on;
plot(1,NovObjPercents,'r*');
plot(2,FamObjPercents,'r*');
plot(3,NovEmptyPercents,'r*');
plot(4,FamEmptyPercents,'r*');
if ischar(minutes)==1
    tempstring='Time Spent in Each Quadrant, All Exposures, All Minutes';
elseif minutes(1)==minutes(2)
    tempstring=['Time Spent in Each Quadrant, All Exposures, Minute ' num2str(minutes(1))];
else
    tempstring=['Time Spent in Each Quadrant, All Exposures, Minutes ' num2str(minutes(1)) ' Through ' num2str(minutes(2))];
end
title(tempstring);
ylabel('Percent of Total Time');
xlabel('Quadrants');
set(gca,'XTickLabel',{'Novel Object' 'Familiar Object' 'Novel Empty' 'Familiar Empty'});
set(gca,'YLim',[0 1]);
hold off;





