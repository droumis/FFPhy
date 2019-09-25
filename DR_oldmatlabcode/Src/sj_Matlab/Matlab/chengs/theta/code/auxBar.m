function auxBar(m, dev, opt, pval)
%function auxBar(m, dev, opt, pval)
%
% Make pretty bar plot with errorbars and comparison brackets
%
% opt.*:  title outname plotsize xlabel ylabel xtick xticklabels xlim ylim 

if nargin<2; dev= []; end
if nargin<3; opt= []; end
if nargin<4; pval= []; end

nx= size(m,1);
ny= size(m,2);

if ~isfield(opt, 'xlabel'); opt.xlabel= ''; end
if ~isfield(opt, 'ylabel'); opt.ylabel= 'y'; end
if ~isfield(opt, 'xticklabels'); opt.xticklabels= num2str([1:nx]'); end
%if ~isfield(opt, 'plotcol'); opt.plotcol= zeros(nx,ny,3); end
if ~isfield(opt, 'title'); opt.title=''; end
if ~isfield(opt, 'outname'); opt.outname=''; end
if ~isfield(opt, 'plotsize'); opt.plotsize='mini'; end

if isempty(dev); dev= zeros(nx,ny); end

maxy= max(max(m+dev));
miny= min(min(m-dev));
pad= (maxy-miny)/30;

figure
h= bar(m);
for j=1:ny
        tmp= get(get(h(j), 'Children'), 'XData');
%        tmp= get(h(j), 'XData');
    xd(:,j)= mean(tmp(2:3,:))';
end
    
cla
hold on
sl= 0.03*maxy;
for i=1:nx
    bw= diff(tmp(2:3,1));
    for j=1:ny
        hs= bar(xd(i,j)+[0 1],[m(i,j) 0], bw);
        if iscell(opt.plotcol)
            set(hs, 'FaceColor', opt.plotcol{j}(i,:));
        else
            set(hs, 'FaceColor', opt.plotcol(i,:));
        end
        errorbar(xd(i,j), m(i,j), dev(i,j), 'Color', 0.4*[1 1 1]);
        if j>1 & ~isempty(pval)
            if pval(i,j-1)<0.05;
    %            yl= max(m(i,j-1:j)+ 1.2*dev(i,j-1:j));
    %            sl= (yl- max(m(i,j-1:j)))*0.9;
                yl= max(m(i,j-1:j)+dev(i,j-1:j)) + 2*sl;
                if(1.1*yl>maxy) maxy= 1.1*yl; end
                hl= line(xd(i,[j-1,j-1,j,j])+[1 1 -1 -1]*bw/8, [yl-sl yl yl yl-sl]);
%                hl= line(xd(i,[j-1,j-1,j,j]), [yl-sl yl yl yl-sl]);
                set(hl, 'Color', 'k');

                if pval(i,j-1)<1e-3;
                    sigstr= '***';
                elseif pval(i,j-1)<1e-2;
                    sigstr= '**';
                elseif pval(i,j-1)<0.05;
                    sigstr= '*';
                else
                    sigstr= '';
                end
                ht= text(mean(xd(i,j-1:j)), yl+sl/2, sigstr);
                set(ht, 'HorizontalAlignment', 'center');
                set(ht, 'VerticalAlignment', 'middle');
            end
        end
    end
end
if isfield(opt, 'ylim') 
    set(gca, 'YLim', opt.ylim);
else
    if miny>0; 
        set(gca, 'YLim', [0 maxy+pad])
    else
        set(gca, 'YLim', [miny-pad maxy+pad])
    end
end
if isfield(opt, 'xlim') 
    set(gca, 'XLim', opt.xlim);
else
    set(gca, 'XLim', [0.5 nx+0.5]);
end
%keyboard
xlabel(opt.xlabel); ylabel(opt.ylabel);
if isfield(opt, 'xtick')
    set(gca, 'XTick', opt.xtick);
else
    set(gca, 'XTick', [1:nx])
end
set(gca, 'XTickLabel', opt.xticklabels);

if ~isempty(opt.outname)
    figname=['bar_' opt.title '_' opt.outname];
else
    figname=['bar_' opt.title];
end
set(gcf, 'Name', figname);
saveas(gcf, figname);
myprint(opt.plotsize, figname);

