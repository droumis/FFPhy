% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel'};

% Epoch selection
epochfilter{1} = 'isequal($type,''run'') & isequal($session,''familiar'')';
epochfilter{2} = 'isequal($type,''run'') & isequal($session,''novel'')';
epochfilter{3} = 'isequal($type,''run'') & isequal($session,''supernovel'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($immobilitytime < 10)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcnovelobjectbehavior', {'task','pos'});
f = runfilter(f);

%% PLOT

% Look at percent time spent in each quadrant type for all minutes
obj = []; nov = []; emp = []; fam = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if f(an).output{d}(e).totalpercenttime(1) ~= -1
                obj = [obj; sum(f(an).output{d}(e).totalpercenttime([1 2]))];
                emp = [emp; sum(f(an).output{d}(e).totalpercenttime([3 4]))];
                nov = [nov; sum(f(an).output{d}(e).totalpercenttime([1 3]))];
                fam = [fam; sum(f(an).output{d}(e).totalpercenttime([2 4]))];
            end
        end
    end
end
invalid = (obj< 0) | (nov < 0) | (emp < 0) | (fam < 0);
obj(invalid)= []; nov(invalid) = []; emp(invalid) = []; fam(invalid) = [];

figure
hold on
bar([1 2 4 5],[mean(obj) mean(emp) mean(nov) mean(fam)],'b')
errorbar2([1 2 4 5],[mean(obj) mean(emp) mean(nov) mean(fam)], ...
    [std(obj)./sqrt(length(obj)-1) std(emp)./sqrt(length(obj)-1) std(nov)./sqrt(length(obj)-1) std(fam)./sqrt(length(obj)-1)],'k')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'Object'},{'Empty'},{'Novel'},{'Familiar'}])
ylabel('Percent time')

% Both obj vs. empty and novel vs. familiar are significant by paired
% ttest and paired signrank for all times
% paired ttest: obj vs. empty p<0.001, nov vs. fam p<0.01
% paired signrank: obj vs. empty p<0.01, nov vs. fam p<0.01 

% Look at percent time spent in each quadrant type for first minute
obj = []; nov = []; emp = []; fam = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if f(an).output{d}(e).totalpercenttime(1) ~= -1
                obj = [obj; sum(f(an).output{d}(e).percenttime([1 2]))];
                emp = [emp; sum(f(an).output{d}(e).percenttime([3 4]))];
                nov = [nov; sum(f(an).output{d}(e).percenttime([1 3]))];
                fam = [fam; sum(f(an).output{d}(e).percenttime([2 4]))];
            end
        end
    end
end
invalid = (obj< 0) | (nov < 0) | (emp < 0) | (fam < 0);
obj(invalid)= []; nov(invalid) = []; emp(invalid) = []; fam(invalid) = [];

figure
hold on
bar([1 2 4 5],[mean(obj) mean(emp) mean(nov) mean(fam)],'b')
errorbar2([1 2 4 5],[mean(obj) mean(emp) mean(nov) mean(fam)], ...
    [std(obj)./sqrt(length(obj)-1) std(emp)./sqrt(length(obj)-1) std(nov)./sqrt(length(obj)-1) std(fam)./sqrt(length(obj)-1)],'k')

set(gca,'xtick',[1 2 4 5],'xticklabel',[{'Object'},{'Empty'},{'Novel'},{'Familiar'}])
ylabel('Percent time')

% Both obj vs. empty and novel vs. familiar for first minute is significant:
% paired ttest: obj vs. empty p<0.001, nov vs. fam p<0.001
% paired signrank: obj vs. empty p<0.001, nov vs. fam p<0.001 

%Look at preference for novelty over time
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
fampref(11:end,:) = []; novpref(6:end,:) = []; suppref(6:end,:) = [];

figure
hold on
plot(1:10,nanmean(fampref,2),'k')
plot(12:16,nanmean(novpref,2),'c')
plot(18:22,nanmean(suppref,2),'m')
legend('Familiarization','Novel Location','Novel Object & Location')
errorbar2(1:10,nanmean(fampref,2),nanstd(fampref,[],2)./sqrt(famcount-1),'k')
errorbar2(12:16,nanmean(novpref,2), nanstd(novpref,[],2)./sqrt(novcount-1),'c')
errorbar2(18:22,nanmean(suppref,2),nanstd(suppref,[],2)./sqrt(supcount-1),'m')
set(gca,'xtick',[1:10 12:16 18:22],'xticklabel',[1:10 1:5 1:5])
xlabel('Cummulative Time (min)')
ylabel('Novelty preference')

% Signrank test for first minute:
% Familiar: p>0.75
% Novel: p<0.01
% SuperNovel: p<0.001

% Signrank test for last minute:
% Familiar: p>0.47
% Novel: p>0.69
% SuperNovel: p<0.001

%Controlling for number of comparisons, p<0.01 to be considered significant
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_preference_over_time.pdf', m, d, y);
print('-dpdf', savestring)



%Look at total preference for novelty
fampref = []; novpref = []; suppref = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if d == 1
                fampref = [fampref f(an).output{d}(e).totalpreference];
            elseif d == 2
                novpref = [novpref f(an).output{d}(e).preference(5)];
            elseif d == 3
                suppref = [suppref f(an).output{d}(e).preference(5)];
            end
        end
    end
end
figure
hold on
bar([1 2 3],[mean(fampref) mean(novpref) mean(suppref)],'b')
errorbar2([1 2 3],[mean(fampref) mean(novpref) mean(suppref)],...
    [std(fampref)./sqrt(length(fampref)-1) std(novpref)./sqrt(length(novpref)-1) std(suppref)./sqrt(length(novpref)-1)],'k')
set(gca,'xtick',[1 2 3],'xticklabel',[{'Familiar'},{'Novel'},{'SuperNovel'}])
ylabel('Novelty preference')