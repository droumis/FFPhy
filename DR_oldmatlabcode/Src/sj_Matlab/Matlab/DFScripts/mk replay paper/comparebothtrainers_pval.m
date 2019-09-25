
%script for comparing the p-values of events from E1 and E2

totalmatrix = [];

for i = 1:size(out,1);
    tmpmatrix = [];
    tmpmatrix = [out{i,1}(:,3:5) out{i,2}(:,3:5)];
    totalmatrix = [totalmatrix; tmpmatrix];
end

totalmatrix = totalmatrix(find((totalmatrix(:,3) >= 5)&(totalmatrix(:,6) >= 5)&(totalmatrix(:,2)>0)),:);
plot(totalmatrix(:,1),totalmatrix(:,4),'.');

numvals = length(totalmatrix);

corrval = corrcoef(totalmatrix(:,1),totalmatrix(:,4));
corrval = corrval(2,1);

for i = 1:10000
    tmpval = [totalmatrix(:,1), totalmatrix(randperm(numvals),4)];
    tmpcorrval = corrcoef(tmpval(:,1),tmpval(:,2));
    permcorrval(i,1) = tmpcorrval(2,1);
end

sum(corrval > permcorrval)/10000