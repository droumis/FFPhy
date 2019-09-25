load '/data13/mcarr/RipplePaper/decodefilterA.mat'
valid_ripples = cell(length(decodefilter),1);
valid_epochs = cell(length(decodefilter),1);
for an = 1:length(decodefilter)
    valid_ripples{an} = cell(length(decodefilter(an).epochs{1}),1);
    valid_epochs{an} = decodefilter(an).epochs{1};
    for d = 1:length(decodefilter(an).epochs{1})
        if ~isempty(decodefilter(an).output{1}(d).eventtime)
            valid_ripples{an}{d}.time = decodefilter(an).output{1}(d).eventtime(~isnan(decodefilter(an).output{1}(d).pvalue),1);
            valid_ripples{an}{d}.pvalue = decodefilter(an).output{1}(d).pvalue(~isnan(decodefilter(an).output{1}(d).pvalue));
        end
    end
end
            
time = f(1).output{1}(1).time;
freq = f(1).output{1}(1).frequency;
c = cell(length(f),1);
for an = 1:length(f)
    animal_count = 0;
    valid_count = 1;
    for d = 1:length(f(an).eegdata{1})
        tmp = [];
        for e = 1:length(f(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if rowfind(f(an).epochs{1}(d,:),valid_epochs{an})>0 && ~isempty(f(an).output{1}(animal_count).ripples)
                if isempty(tmp)
                    ind = lookup(valid_ripples{an}{valid_count}.time,f(an).output{1}(animal_count).ripples);
                    tmp = f(an).output{1}(animal_count).coherence(:,:,ind);
                    tet_count = 1;
                else
                    tmp = tmp + f(an).output{1}(animal_count).coherence(:,:,ind);
                    tet_count = tet_count+1;
                end
            end       
        end
        c{an}{valid_count} = tmp./tet_count;
        if rowfind(f(an).epochs{1}(d,:),valid_epochs{an})>0
            valid_count = valid_count+1;
        end
    end
end

%Look at gamma coherence for significant and nonsignificant replay events
subs = [0.05 1];
gam = []; p = [];
for an = 1%:length(c)
    for d = 1:length(c{an})
        if isempty(gam) && ~isempty(c{an}{d})
            gam = c{an}{d};
            p = valid_ripples{an}{d}.pvalue;
        elseif ~isempty(c{an}{d})
            gam = cat(3,gam,c{an}{d});
            p = [p; valid_ripples{an}{d}.pvalue];
        end
    end
end
p =lookup(p,subs,1);
a = zeros(length(subs),length(time));   v = zeros(length(subs),length(time));

for i = 1:size(gam,1)
    a(:,i) = accumarray(p,squeeze(mean(gam(i,2:9,:),2)),[length(subs) 1],@(x) mean(x));
    v(:,i) = accumarray(p,squeeze(mean(gam(i,2:9,:),2)),[length(subs) 1],@(x) std(x)./sqrt(length(x)-1));
end
figure
fill([time time(end:-1:1)],[a(1,:)+v(1,:) a(1,end:-1:1)-v(1,end:-1:1)],'r','Edgecolor','None')
hold on
fill([time time(end:-1:1)],[a(end,:)+v(end,:) a(end,end:-1:1)-v(end,end:-1:1)],'k','Edgecolor','None')
xlabel('Time since CA1 ripple (sec)')
ylabel('CA1-CA3 Gamma Coherence')


%Look at correlation between replay significance and coherence
subs = [0 0.0001 0.01 0.05 0.1 0.5 1];
gam = []; p = [];
for an = 1:length(c)
    for d = 1:length(c{an})
        if isempty(gam) && ~isempty(c{an}{d})
            gam = c{an}{d};
            p = valid_ripples{an}{d}.pvalue;
        elseif ~isempty(c{an}{d})
            gam = cat(3,gam,c{an}{d});
            p = [p; valid_ripples{an}{d}.pvalue];
        end
    end
end
p = lookup(p,subs,1);
a = zeros(length(subs),length(time));
for i = 1:size(gam,1)
    a(:,i) = accumarray(p,squeeze(mean(gam(i,2:9,:),2)),[length(subs) 1],@(x) mean(x));
end
figure
imagesc(time,[],a)
set(gca,'clim',[.4 0.75],'ytick',1:length(subs),'yticklabel',subs)
ylabel('pvalue')
xlabel('Time (sec)')
ylabel('Replay p-value')
colorbar
title('CA1-CA3 gamma coherence')


[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA1CA3_gamma_coherence_pvalue.pdf',m,d,y);
print('-dpdf', savestring)