%function showLinCorr()
%%%key showLinCorr

plotmap=[1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];
binwidth= .2;

passes=1:16;
load (fmaux.select)
meanr= zeros(max(passes),1);
stdr= zeros(max(passes),1);
nr= zeros(max(passes),1);
binloc=[-1+binwidth/2:binwidth:1-binwidth/2];
figure(1);
for i= passes
    s{i}= collectStats('LinCorr',i,4);
%    s{i}.val

%    for n=1:length(s{i}.val)
%        selid= getSelectId(s{i}.cellnum(n,:));
%        traj= select.a{selid}{s{i}.aindex(n)}.traj;
%        fprintf(1, '[%d %d %d %d], traj= %d, r= %f\n', s{i}.cellnum(n,:), ...
%        traj, s{i}.val(n));
%    end
    if ~isfield(s{i},'val'); continue; end
    h=hist(s{i}.val, binloc);
    subplot(4,4, plotmap(i));
    bar(binloc,h);
    hh = findobj(gca,'Type','patch');
    set(hh,'FaceColor','k','EdgeColor','w')
    axis([-1 1 0 max(h)]);
    title(['pass ' num2str(i)]);
    meanr(i)= mean(s{i}.val);
    nr(i)= length(s{i}.val);
    stdr(i)= std(s{i}.val);
    fprintf(1, 'pass %d: mean= %f, std= %f, n= %d\n', i, mean(s{i}.val), std(s{i}.val), length(s{i}.val)); 
end

figure(2)
errorbar(passes, meanr(passes), stdr(passes)./nr(passes));
title([fmaux.prefix]);
axis tight
xlabel('pass');
ylabel('lin corr coeff');

printout= 0;
if printout
    for f=1:2
        figure(f);
        orient landscape;
        print -dps -P811
    end
end
