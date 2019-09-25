%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
epochfilter{2} = 'isequal($type,''sleep'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcconcurrentrip', {'ripples','cellinfo'});
f = runfilter(f);

%save('/data13/mcarr/RipplePaper/concurrentripples.mat','f')


thresh3 = []; thresh5 = []; thresh7 = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            thresh3 = [thresh3; f(an).output{d}(e).concurrent(1)];
            thresh5 = [thresh5; f(an).output{d}(e).concurrent(2)];
            thresh7 = [thresh7; f(an).output{d}(e).concurrent(3)];
        end
    end
end
invalid = isnan(thresh3); thresh3(invalid) = [];
invalid = isnan(thresh5); thresh5(invalid) = [];
invalid = isnan(thresh7); thresh7(invalid) = [];


figure
bar([1 2 3],[nanmean(thresh3) nanmean(thresh5) nanmean(thresh7)],1,'b')
hold on
errorbar2([1 2 3],[nanmean(thresh3) nanmean(thresh5) nanmean(thresh7)],[nanstd(thresh3)./sqrt(sum(~isnan(thresh3))-1) nanstd(thresh5)./sqrt(sum(~isnan(thresh5))-1) nanstd(thresh7)./sqrt(sum(~isnan(thresh7))-1)],'k')
set(gca,'xtick',[1 2 3],'xlim',[0.25 3.75],'xticklabel',[{'>3std'}, {'>5std'}, {'>7std'}])
ylabel('Proportion of concurrent SWRs')
xlabel('SWR amplitude')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_concurrent_SWRs_across_hemispheres.pdf', m, d, y);
print('-dpdf', savestring)
