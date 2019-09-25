% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel'};

% Epoch selection
epochfilter{1} = 'isequal($type,''run'') & isequal($session,''familiar'')';
epochfilter{2} = 'isequal($type,''run'') & isequal($session,''novel'')';
epochfilter{3} = 'isequal($type,''run'') & isequal($session,''supernovel'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($immobilitytime < 5)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcnovelobjectbehavior', {'task','pos'});
f = runfilter(f);

%% Look at percent time spent in object and empty quadrants
objnov = []; empnov = []; objsup = []; empsup = []; objfam = []; empfam = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if d==1
                objfam = [objfam; f(an).output{d}(e).percenttime(2) f(an).output{d}(e).totalpercenttime(2)];
                empfam = [empfam; f(an).output{d}(e).percenttime(4) f(an).output{d}(e).totalpercenttime(4)];
            elseif d==2
                objnov = [objnov; sum(f(an).output{d}(e).percenttime([1 2])) sum(f(an).output{d}(e).totalpercenttime([1 2]))];
                empnov = [empnov; sum(f(an).output{d}(e).percenttime([3 4])) sum(f(an).output{d}(e).totalpercenttime([3 4]))];
            elseif d==3
                objsup = [objsup; sum(f(an).output{d}(e).percenttime([1 2])) sum(f(an).output{d}(e).totalpercenttime([1 2]))];
                empsup = [empsup; sum(f(an).output{d}(e).percenttime([3 4])) sum(f(an).output{d}(e).totalpercenttime([3 4]))];
            end
        end
    end
end

invalid = (objsup(:,1) < 0) | (empsup(:,1) < 0); objsup(invalid,:)= []; empsup(invalid,:) = [];
objsup(objsup == 0) = 0.01; empsup(empsup==0) = 0.01;


%All time
means = [mean(objfam(:,2)) mean(empfam(:,2)) mean(objnov(:,2)) mean(empnov(:,2)) mean(objsup(:,2)) mean(empsup(:,2))];
errors = [std(objfam(:,2)) std(empfam(:,2)) std(objnov(:,2)) std(empnov(:,2)) std(objsup(:,2)) std(empsup(:,2))];
errors = errors./sqrt([size(objfam,1)-1 size(empfam,1)-1 size(objnov,1)-1 size(empnov,1)-1 size(objsup,1)-1 size(empsup,1)-1]);
figure
hold on
bar([1 2],means([1 2]),'b')
bar([4 5],means([3 4]),'r')
bar([7 8],means([5 6]),'c')
legend('Familiar','Novel','SuperNovel')
errorbar2([1 2 4 5 7 8],means,errors,'k')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',[{'Object'},{'Empty'},{'Object'},{'Empty'},{'Object'},{'Empty'}])
ylabel('Percent time spent in quadrant')
title('All minutes')

p = signrank(objfam(:,2)./empfam(:,2),1); %pvalue < 0.001
p = signrank(objnov(:,2)./empnov(:,2),1); %pvalue < 0.05
p = signrank(objsup(:,2)./empsup(:,2),1); %pvalue > 0.25

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_object_vs_empty_allminutes.pdf', m, d, y);
print('-dpdf', savestring)


%First minute
means = [mean(objfam(:,1)) mean(empfam(:,1)) mean(objnov(:,1)) mean(empnov(:,1)) mean(objsup(:,1)) mean(empsup(:,1))];
errors = [std(objfam(:,1)) std(empfam(:,1)) std(objnov(:,1)) std(empnov(:,1)) std(objsup(:,1)) std(empsup(:,1))];
errors = errors./sqrt([size(objfam,1)-1 size(empfam,1)-1 size(objnov,1)-1 size(empnov,1)-1 size(objsup,1)-1 size(empsup,1)-1]);
figure
hold on
bar([1 2],means([1 2]),'b')
bar([4 5],means([3 4]),'r')
bar([7 8],means([5 6]),'c')
legend('Familiar','Novel','SuperNovel')
errorbar2([1 2 4 5 7 8],means,errors,'k')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',[{'Object'},{'Empty'},{'Object'},{'Empty'},{'Object'},{'Empty'}])
ylabel('Percent time spent in quadrant')
title('First minute')

p = signrank(objfam(:,1)./empfam(:,1),1); %pvalue < 1e-5
p = signrank(objnov(:,1)./empnov(:,1),1); %pvalue < 0.001
p = signrank(objsup(:,1)./empsup(:,1),1); %pvalue > 0.001

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_object_vs_empty_firstminute.pdf', m, d, y);
print('-dpdf', savestring)


%% Look at preference for novelty over time
fampref = NaN(30,100); novpref = NaN(30,100); suppref = NaN(30,100);
famcount = 1; novcount = 1; supcount = 1;
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if d == 1
                fampref(1:length(f(an).output{d}(e).preference),famcount) = f(an).output{d}(e).preference;
                famcount = famcount+1;
            elseif d == 2
                novpref(1:length(f(an).output{d}(e).preference),novcount) = f(an).output{d}(e).preference;
                novcount = novcount+1;
            elseif d == 3
                suppref(1:length(f(an).output{d}(e).preference),supcount) = f(an).output{d}(e).preference;
                supcount = supcount+1;
            end
        end
    end
end
fampref(:,famcount:end) = []; novpref(:,novcount:end) = [];  suppref(:,supcount:end) = [];
fampref(6:end,:) = []; novpref(6:end,:) = []; suppref(6:end,:) = [];

figure
hold on
plot(1:5,nanmean(fampref,2),'k',1:5,nanmean(novpref,2),'c',1:5,nanmean(suppref,2),'m')
legend('Familiarization','Novel Location','Novel Object & Location')
errorbar2(1:5,nanmean(fampref,2),nanstd(fampref,[],2)./sqrt(famcount-1),'k')
errorbar2(1:5,nanmean(novpref,2), nanstd(novpref,[],2)./sqrt(novcount-1),'c')
errorbar2(1:5,nanmean(suppref,2),nanstd(suppref,[],2)./sqrt(supcount-1),'m')
set(gca,'xtick',1:10)
xlabel('Cummulative Time (min)')
ylabel('Novelty preference')

% Rank sum test for first minute:
% Familiar vs. Novel: p<0.05
% Familiar vs. SuperNovel: p<0.01

% Rank sum test for last minute:
% Familiar vs. Novel: p>0.85
% Familiar vs. SuperNovel: p<0.01

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_preference_over_time.pdf', m, d, y);
print('-dpdf', savestring)


%Look at total preference for novelty
tmp1 = fampref(~isnan(fampref(1:5,:)));
tmp2 = novpref(~isnan(novpref(1:5,:)));
tmp3 = suppref(~isnan(suppref(1:5,:)));
means = [mean(tmp1) mean(tmp2) mean(tmp3)];
errors = [std(tmp1) std(tmp2) std(tmp3)];
errors = errors./sqrt([length(tmp1)-1 length(tmp2)-1 length(tmp3)-1]);

figure
hold on
bar(1,means(1),'b')
bar(2,means(2),'r')
bar(3,means(3),'c')
errorbar2([1 2 3],means,errors,'k')
set(gca,'xtick',[1 2 3],'xticklabel',[{'Familiar'},{'Novel'},{'SuperNovel'}])
ylabel('Novelty preference')

p = signrank(tmp1); %pvalue > 0.3
p = signrank(tmp2); %pvalue < 0.05
p = signrank(tmp3); %pvalue < 1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_preference_alltime.pdf', m, d, y);
print('-dpdf', savestring)


