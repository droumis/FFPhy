%% plots num cells, day by day (from first epoch)

%% ingredient: tetinfo


figure
hold on
e=1;   % pick arbitrary epoch with a day

for d=3:length(tetinfo)
    dummy_numcells=nan(1,length(tetinfo{d}{e}));
    dummy_area=cell(1,length(tetinfo{d}{e}));
    for t=1:length(tetinfo{d}{e})
        if ~isempty(tetinfo{d}{e}{t})
            dummy_numcells(t)=tetinfo{d}{e}{t}.numcells;
            dummy_area{t}=tetinfo{d}{e}{t}.area;
        else
            dummy_numcells(t)=0;
            dummy_area{t}='none';
        end
    end
    subplot(length(tetinfo),1,d)
    bar(dummy_numcells);
    set(gca,'XTick',1:25,'XTickLabel',dummy_area);
end