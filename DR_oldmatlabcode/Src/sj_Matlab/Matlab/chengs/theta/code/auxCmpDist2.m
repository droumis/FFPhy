function auxCmpDist2(x, sel, opt)
%function auxCmpDist2(x, sel, opt)
%
% compare samples in x{i}
%
%   sel{i}     description of samples
%       label
%
%   opt.plot='cdf'; {'cdf', 'pdf', 'hist', 'bar'}
%   opt.name        name of variable for internal use
%   opt.title       name of variable w/o special char
%   opt.label       name of variable for labeling in figure
%   opt.ref= 1
%   opt.plotcol= {'k', 'r', ...}
%   opt.size
%   opt.nbins

if nargin<3; opt= []; end
if ~isfield(opt, 'title'); opt.title='comparison'; end
if ~isfield(opt, 'ref'); opt.ref=1; end
if ~isfield(opt, 'plot'); opt.plot='cdf'; end
if ~isfield(opt, 'plotcol');
    plotcol(1:6,:)=[[0 0 0]; [1 0 0]; [0 0 1]; [0 1 0]; [.7 0 1]; [0 1 1]];
else
    plotcol= opt.plotcol;
end
nbins= 10;      if isfield(opt, 'nbins') nbins= opt.nbins; end
norm= 0;        if isfield(opt, 'norm');  norm= opt.norm; end
plotsize= 'small'; if isfield(opt, 'size');  plotsize= opt.size; end
showlegend= 0;  if isfield(opt, 'legend');  showlegend= opt.legend; end
outname= '';  if isfield(opt, 'outname');  outname=  opt.outname; end

nsel= length(x);

fprintf(1, '\ncomparing %s %s', opt.label, outname);
fprintf(1, '\t(ranksum test, KS test, t-test)\n');
for is=1:nsel
    if is==opt.ref; continue; end
    if length(x{opt.ref})<2 | length(x{is})<2
        prank= nan;
        pKS= nan;
        pT= nan;
    else
        prank= ranksum(x{opt.ref}, x{is});
        [h,pKS]= kstest2(x{opt.ref}, x{is});
        [h,pT]= ttest2(x{opt.ref}, x{is});
    end
    fprintf(1, '  "%s" (n=%d) vs. "%s" (n=%d)\t(p=%.2g, p=%.2g, p=%.2g)\n', ...
        sel{opt.ref}.label, length(x{opt.ref}), sel{is}.label, length(x{is}),...
        prank, pKS, pT);
end

for is=1:nsel
    fprintf(1, '%s %s: %s= %.2f +/- %.2f\n', outname, sel{is}.label, ...
        opt.label, nanmean(x{is}), nanstd(x{is}));
end


if strcmp(opt.plot, 'none'); return; end

%figure;
hold on
lstr= {};
nsel= length(x);

minv= nan; maxv= nan;
for is=1:nsel
    if size(x{is},2)>1; x{is}= x{is}'; end
    minv= min([minv; x{is}]);
    maxv= max([maxv; x{is}]);
end

showDist= 1;

switch opt.plot
case {'hist', 'pdf'}

    if isfield(opt, 'xrange')
        borders= linspace(opt.xrange(1), opt.xrange(2), nbins+1);
    else
        borders= linspace(.95*minv, 1.05*maxv, nbins+1);
    end
    pb= borders(1:end-1)+(borders(2)-borders(1))/2;

    for is=1:nsel
        n= histc(x{is},borders);
        if norm; n= n/sum(n); end
        Nis(:,is)= n(1:end-1);
    end
    switch opt.plot
    case 'hist'
        h= bar(pb,Nis);
        for is=1:nsel; set(h(is), 'FaceColor', plotcol(is,:)); end
    case 'pdf'
        h= plot(pb,Nis);
        for is=1:nsel; set(h(is), 'Color', plotcol(is,:)); end
    end

    if norm
        ylabel('fraction');
    else
        ylabel('count');
    end
    maxy= max(max(Nis));

    if isfield(opt, 'yrange')
        set(gca, 'YLim', opt.yrange);
    else
        set(gca, 'YLim', [0, 1.05*maxy]);
    end
case {'cdf', 'cdf-p'}
    for is=1:nsel
        h= cdfplot(x{is});
        set(h, 'Color', plotcol(is,:));
%        set(h, 'LineStyle', '-.');
        hold on
    end
    title('');
    ylabel('cum fraction');
case 'bar'
    showDist= 0;
    for i=1:length(x)
        Zm(i,1)= nanmean(x{i});
        Zdev(i,1)= nanstd(x{i})/sqrt(sum(isfinite(x{i})));
    end
    for is=1:nsel
        lstr{is}= sel{is}.label;
    end
    opt.xticklabels= lstr;
    figure
    auxBar(Zm, Zdev, opt);
end

%keyboard
if showDist
    if isfield(opt, 'xrange')
        xrange= opt.xrange;
    else
        xrange= [0.95*minv, 1.05*maxv];
    end
    set(gca, 'XLim', xrange);

    if isfield(opt, 'label'); xlabel(opt.label); end
else
    if isfield(opt, 'label'); ylabel(opt.label); end
end


if strfind(opt.plot, '-p') & nsel==2
    tstr= sprintf('p=%.2g', prank);
    ht= text(xrange(2), 0, tstr);
    set(ht,'VerticalAlignment', 'bottom')
    set(ht,'HorizontalAlignment', 'right')
end

if showlegend 
    for is=1:nsel
        lstr{is}= sel{is}.label;
    end
    legend(lstr, 0);
end

%if strfind(var, 'PeakRatio') set(gca, 'xlim', [0,10]); end
hold off
if ~isempty(outname)
    figname=[opt.plot '_' opt.title '_' outname];
else
    figname=[opt.plot '_' opt.title];
end
set(gcf, 'Name', figname);
myprint(plotsize, figname);
