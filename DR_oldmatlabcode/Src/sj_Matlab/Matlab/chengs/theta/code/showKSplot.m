function [violations, viol_ind]= showKSplot(num, direct, mode)
%function [violations, viol_ind]= showKSplot(num, direct, mode)
%  Show KS plot for rescaled isi
%  
%  direct=1     1: Show direct model estimates. 0: show stats calculation
%  mode= 0      how to collect and display results
%                >0: add results to persistent variables, do not display
%                 0: analyze, display and print results
%                -1: display and print results, do not analyze

global KSdist color makegraph fmaux

color=     1;
prettyPlot= 1;


if nargin < 2 | isempty(direct); direct= 1; end
if nargin < 3 | isempty(mode); mode=0; end

global violations viol_ind
if mode==0 | mode== 1
    violations= []; viol_ind= [];
end

if mode== 0;
    makegraph= 1;
else 
    makegraph= 0;
end

if mode >=0
    if nargin < 1 | isempty(num)
        [d,e,t,c,ncells]= startCellList;
        viol= [0 0];  % number of KS test violation at the 95%- and 99% level
        vind= [];     % whether cells that violate KS test at 95% and 99% level
        n= 0;
        while ~isempty(d)
            v= auxrun(d,e,t,c, direct);
            viol= viol+v;
            n= n+1;
            viol_ind(end+1,1)= fmaux.currentCell; 
            violations(end+1,:)= v;
            [d,e,t,c]= getNextCell;
        end
        fprintf(1, '%d (%d) out of %d (%.1f%%/ %.1f%%) within 95%% (99%%)-CI.\n', n-viol, n, (n-viol)/n*100);
    else
        d=num(1); e=num(2); t=num(3); c=num(4);
        v= auxrun(d,e,t,c, direct);
        vind= v;
    end
end

if mode <= 0
    n= size(violations,1);
    nviol= sum(violations);
    fprintf(1, '%d (%d) out of %d (%.1f%%/ %.1f%%) within 95%% (99%%)-CI.\n',...
    n-nviol, n, 100*(n-nviol)/n);

end


function viol= auxrun(d,e,t,c, direct)

KSdist95 = 1.36;
KSdist99 = 1.63;
global sorted

if direct
    global KSdist fmaux prettyPlot color adaptest makegraph
    loadVar('.','adaptest',d);
    sorted= sort(adaptest{d}{e}{t}{c}.model.rescaled_isi);
else
    global KSdist fmaux prettyPlot color stats makegraph
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

if (sum( abs([1:N]'/N-sorted) > KSdist95/sqrt(N)))
    viol(1)= 1;
    fprintf(1, '  outside 95%%-CI. \t');
else 
    viol(1)= 0;
    fprintf(1, '  ok! inside 95%%-CI. \t');
end
if (sum( abs([1:N]'/N-sorted) > KSdist99/sqrt(N)))
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
    if ~prettyPlot
        subplot(2,1,1);
        %plot(b, h);
        hist(sorted,20)
        hh = findobj(gca,'Type','patch');
        set(hh,'FaceColor','k','EdgeColor','w')
        subplot(2,1,2);
    end

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
    if ~prettyPlot
        tstring= sprintf('d= %d, e= %d, t= %d, c= %d. viol= [%d %d]', d,e,t,c, viol);
        title(tstring);
        orient tall
    end
    figname= sprintf('KSplot-%d-%d-%d-%d', d,e,t,c);
    xlabel('cum. frac.')
    ylabel('cum. frac.')
    set(gcf, 'Name', figname);
    myprint('small', figname);
%    saveas(gcf, figname);
end

