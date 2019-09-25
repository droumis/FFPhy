function out=multiepochripthresh(tetrode,minutes)
% combines the ripple threshold per quad data across multiple epochs, correcting for the
% rotation of objects.  
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


% get ripple threshold by quad data for each desired epoch
[threshsA, typesA]=getripthresh(dayA,epochA,tetrode,minutes);
[threshsB, typesB]=getripthresh(dayB,epochB,tetrode,minutes);
[threshsC, typesC]=getripthresh(dayC,epochC,tetrode,minutes);
[threshsD, typesD]=getripthresh(dayD,epochD,tetrode,minutes);
close all; %to close the figures for each epoch, can always comment this out if want to see them


% find which quadrant was the NovObj quad in each epoch
whereNOa=strcmp('NovObj',typesA); %where is the NovObj in epoch A
whereNOb=strcmp('NovObj',typesB);
whereNOc=strcmp('NovObj',typesC);
whereNOd=strcmp('NovObj',typesD);
% find the ripple threshold in NovObj quad for each epoch and put those
% together. 
NovObjThresh=[0 0 0 0];
NovObjThresh(1)=threshsA(1,whereNOa);
NovObjThresh(2)=threshsB(1,whereNOb);
NovObjThresh(3)=threshsC(1,whereNOc);
NovObjThresh(4)=threshsD(1,whereNOd);
NovObjAvg=nanmean(NovObjThresh);

% find which quadrant was the FamObj quad in each epoch
whereFOa=strcmp('FamObj',typesA); %where is the FamObj in epoch A
whereFOb=strcmp('FamObj',typesB);
whereFOc=strcmp('FamObj',typesC);
whereFOd=strcmp('FamObj',typesD);
% find the ripple threshold in FamObj quad for each epoch and put those
% together. 
FamObjThresh=[0 0 0 0];
FamObjThresh(1)=threshsA(1,whereFOa);
FamObjThresh(2)=threshsB(1,whereFOb);
FamObjThresh(3)=threshsC(1,whereFOc);
FamObjThresh(4)=threshsD(1,whereFOd);
FamObjAvg=nanmean(FamObjThresh);

% find which quadrant was the NovEmpty quad in each epoch
whereNEa=strcmp('NovEmpty',typesA); %where is the NovEmpty in epoch A
whereNEb=strcmp('NovEmpty',typesB);
whereNEc=strcmp('NovEmpty',typesC);
whereNEd=strcmp('NovEmpty',typesD);
% find the ripple threshold in NovEmpty quad for each epoch and put those
% together. 
NovEmptyThresh=[0 0 0 0];
NovEmptyThresh(1)=threshsA(1,whereNEa);
NovEmptyThresh(2)=threshsB(1,whereNEb);
NovEmptyThresh(3)=threshsC(1,whereNEc);
NovEmptyThresh(4)=threshsD(1,whereNEd);
NovEmptyAvg=nanmean(NovEmptyThresh);

% find which quadrant was the FamEmpty quad in each epoch
whereFEa=strcmp('FamEmpty',typesA); %where is the FamEmpty in epoch A
whereFEb=strcmp('FamEmpty',typesB);
whereFEc=strcmp('FamEmpty',typesC);
whereFEd=strcmp('FamEmpty',typesD);
% find the ripple threshold in FamEmpty quad for each epoch and put those
% together. 
FamEmptyThresh=[0 0 0 0];
FamEmptyThresh(1)=threshsA(1,whereFEa);
FamEmptyThresh(2)=threshsB(1,whereFEb);
FamEmptyThresh(3)=threshsC(1,whereFEc);
FamEmptyThresh(4)=threshsD(1,whereFEd);
FamEmptyAvg=nanmean(FamEmptyThresh);


out=[NovObjAvg FamObjAvg NovEmptyAvg FamEmptyAvg];

figure;
h=bar(out);
hold on;
plot(1,NovObjThresh,'r*');
plot(2,FamObjThresh,'r*');
plot(3,NovEmptyThresh,'r*');
plot(4,FamEmptyThresh,'r*');
if ischar(minutes)==1
    tempstring='Average Threshold of Ripples in Each Quadrant, All Exposures, All Minutes';
elseif minutes(1)==minutes(2)
    tempstring=['Average Threshold of Ripples in Each Quadrant, All Exposures, Minute ' num2str(minutes(1))];
else
    tempstring=['Average Threshold of Ripples in Each Quadrant, All Exposures, Minutes ' num2str(minutes(1)) ' Through ' num2str(minutes(2))];
end
title(tempstring);
ylabel('Threshold (in std)');
xlabel('Quadrants');
set(gca,'XTickLabel',{'Novel Object' 'Familiar Object' 'Novel Empty' 'Familiar Empty'});
hold off;
