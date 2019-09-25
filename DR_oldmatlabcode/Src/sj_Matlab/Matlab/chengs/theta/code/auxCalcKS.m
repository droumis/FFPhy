function [viol,KS]= auxCalcKS(d,e,t,c, direct, makegraph)

KSdist95 = 1.36;
KSdist99 = 1.63;
color= 1;

if direct
    global KSdist fmaux adaptest 
    loadVar('.','adaptest',d);
    sorted= sort(adaptest{d}{e}{t}{c}.model.rescaled_isi);
else
    global KSdist fmaux stats makegraph
    loadFile(['stats-' fmaux.selectid '.mat']);
    sorted= sort(stats{d}{e}{t}{c}.RescaledSpikes{1});
end
neg= find(sorted < 0);
if length(neg)/ length(sorted) > .01
    warning(sprintf('too many invalid rescaled isi''s %d/%d', ...
        length(neg), length(sorted)));
end
sorted(neg)= 0;

N= length(sorted);
KS= max(abs([1:N]'/N-sorted));

if KS> KSdist95/sqrt(N)
    viol(1)= 1;
    fprintf(1, '  outside 95%%-CI. \t');
else 
    viol(1)= 0;
    fprintf(1, '  ok! inside 95%%-CI. \t');
end
if KS> KSdist99/sqrt(N)
    viol(2)= 1;
    fprintf(1, '  outside 99%%-CI.\n');
else 
    viol(2)= 0;
    fprintf(1, '  ok! inside 99%%-CI.\n');
end
%    using Matlab's K-S test give same result (cross-check)
%    viol= kstest(sorted,([0:N]'/N)*[1 1]) 

if(makegraph)
    figure
    subplot(2,1,1);
    %plot(b, h);
    hist(sorted,20)
    hh = findobj(gca,'Type','patch');
    set(hh,'FaceColor','k','EdgeColor','w')
    subplot(2,1,2);

    if (color)
        % Find KS statistic in uniform domain
        plot(sorted, ([1:N])/N, 'b', 'LineWidth', 2);  
        hold on
        plot(0:.01:1,0:.01:1, 'g', 'LineWidth', 2); 
        %    plot(0:.01:1, [0:.01:1]+KSdist/sqrt(N), 'r', 0:.01:1,[0:.01:1]-KSdist/sqrt(N), 'r' ); 
        plot(0:.01:1, [0:.01:1]+KSdist95/sqrt(N), 'r', 0:.01:1,[0:.01:1]-KSdist95/sqrt(N), 'r' ); 
        plot(0:.01:1, [0:.01:1]+KSdist99/sqrt(N), 'r', 0:.01:1,[0:.01:1]-KSdist99/sqrt(N), 'r' ); 
        hold off
    else
        plot(sorted, ([1:N])/N, 'k', 'LineWidth', 2.0);  
        hold on
        h = plot(0:.01:1,0:.01:1, 'LineWidth', 1); 
        set(h, 'Color', [.4 .4 .4]);
        h = plot(0:.01:1, [0:.01:1]+KSdist/sqrt(N), '--', 0:.01:1,[0:.01:1]-KSdist/sqrt(N), '--' ); 
        set(h(1), 'Color', [.4 .4 .4]);
        set(h(2), 'Color', [.4 .4 .4]);
        hold off
    end
    axis([0 1 0 1]); axis square
    tstring= sprintf('d= %d, e= %d, t= %d, c= %d. viol= [%d %d]', d,e,t,c, viol);
    title(tstring);
    orient tall
    figname= sprintf('KSplot-%d-%d-%d-%d', d,e,t,c);
%    print('-dpng', figname);
%    saveas(gcf, figname);
end

