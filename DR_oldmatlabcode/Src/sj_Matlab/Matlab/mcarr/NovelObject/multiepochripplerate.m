function out=multiepochripplerate(minutes)
% combines the ripple rate data across multiple epochs, correcting for the
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


% get ripple rate data for each desired epoch
[ripratesA, typesA]=getripplerate(dayA,epochA,minutes);
[ripratesB, typesB]=getripplerate(dayB,epochB,minutes);
[ripratesC, typesC]=getripplerate(dayC,epochC,minutes);
[ripratesD, typesD]=getripplerate(dayD,epochD,minutes);
close all; %to close the figures for each epoch, can always comment this out if want to see them


% find which quadrant was the NovObj quad in each epoch
whereNOa=strcmp('NovObj',typesA); %where is the NovObj in epoch A
whereNOb=strcmp('NovObj',typesB);
whereNOc=strcmp('NovObj',typesC);
whereNOd=strcmp('NovObj',typesD);
% find the ripplerate in NovObj quad for each epoch and put those
% together. 
NovObjRipRate=[0 0 0 0];
NovObjRipRate(1)=ripratesA(1,whereNOa);
NovObjRipRate(2)=ripratesB(1,whereNOb);
NovObjRipRate(3)=ripratesC(1,whereNOc);
NovObjRipRate(4)=ripratesD(1,whereNOd);
NovObjAvg=nanmean(NovObjRipRate);

% find which quadrant was the FamObj quad in each epoch
whereFOa=strcmp('FamObj',typesA); %where is the FamObj in epoch A
whereFOb=strcmp('FamObj',typesB);
whereFOc=strcmp('FamObj',typesC);
whereFOd=strcmp('FamObj',typesD);
% find the ripplerate in FamObj quad for each epoch and put those
% together. 
FamObjRipRate=[0 0 0 0];
FamObjRipRate(1)=ripratesA(1,whereFOa);
FamObjRipRate(2)=ripratesB(1,whereFOb);
FamObjRipRate(3)=ripratesC(1,whereFOc);
FamObjRipRate(4)=ripratesD(1,whereFOd);
FamObjAvg=nanmean(FamObjRipRate);

% find which quadrant was the NovEmpty quad in each epoch
whereNEa=strcmp('NovEmpty',typesA); %where is the NovEmpty in epoch A
whereNEb=strcmp('NovEmpty',typesB);
whereNEc=strcmp('NovEmpty',typesC);
whereNEd=strcmp('NovEmpty',typesD);
% find the ripplerate in NovEmpty quad for each epoch and put those
% together. 
NovEmptyRipRate=[0 0 0 0];
NovEmptyRipRate(1)=ripratesA(1,whereNEa);
NovEmptyRipRate(2)=ripratesB(1,whereNEb);
NovEmptyRipRate(3)=ripratesC(1,whereNEc);
NovEmptyRipRate(4)=ripratesD(1,whereNEd);
NovEmptyAvg=nanmean(NovEmptyRipRate);

% find which quadrant was the FamEmpty quad in each epoch
whereFEa=strcmp('FamEmpty',typesA); %where is the FamEmpty in epoch A
whereFEb=strcmp('FamEmpty',typesB);
whereFEc=strcmp('FamEmpty',typesC);
whereFEd=strcmp('FamEmpty',typesD);
% find the ripplerate in FamEmpty quad for each epoch and put those
% together. 
FamEmptyRipRate=[0 0 0 0];
FamEmptyRipRate(1)=ripratesA(1,whereFEa);
FamEmptyRipRate(2)=ripratesB(1,whereFEb);
FamEmptyRipRate(3)=ripratesC(1,whereFEc);
FamEmptyRipRate(4)=ripratesD(1,whereFEd);
FamEmptyAvg=nanmean(FamEmptyRipRate);


out=[NovObjAvg FamObjAvg NovEmptyAvg FamEmptyAvg];

figure;
h=bar(out);
hold on;
plot(1,NovObjRipRate,'r*');
plot(2,FamObjRipRate,'r*');
plot(3,NovEmptyRipRate,'r*');
plot(4,FamEmptyRipRate,'r*');
if ischar(minutes)==1
    tempstring='Ripple Rate in Each Quadrant, All Exposures, All Minutes';
elseif minutes(1)==minutes(2)
    tempstring=['Ripple Rate in Each Quadrant, All Exposures, Minute ' num2str(minutes(1))];
else
    tempstring=['Ripple Rate in Each Quadrant, All Exposures, Minutes ' num2str(minutes(1)) ' Through ' num2str(minutes(2))];
end
title(tempstring);
ylabel('Ripples/Second');
xlabel('Quadrants');
set(gca,'XTickLabel',{'Novel Object' 'Familiar Object' 'Novel Empty' 'Familiar Empty'});
hold off;






