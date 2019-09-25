function auxShowStats(vars, opt, analist)
%function auxShowStats(vars, opt, analist)

doprint=0; dopro_retro=0; 

nVars= length(vars);
nOpt= length(opt);

for iO= 1:nOpt
    if isfield(opt{iO}, 'print') & opt{iO}.print; doprint=1 ; end
    if isfield(opt{iO}, 'pro_retro') & opt{iO}.pro_retro; dopro_retro=1 ; end

    if ~isfield(opt{iO}, 'action') error('field opt.action must be set'); end
    switch opt{iO}.action
    case 'compare'
        if dopro_retro
            auxPlotCmp(vars, opt{iO}, 'pro', analist.proind)
            auxPlotCmp(vars, opt{iO}, 'retro', analist.retroind)
        else
            auxPlotCmp(vars, opt{iO}, 'all')
        end
    case 'hist'
        if dopro_retro
            auxHist(vars, opt{iO}, 'pro', analist.proind)
            auxHist(vars, opt{iO}, 'retro', analist.retroind)
        else
            auxHist(vars, opt{iO}, 'all')
        end
    otherwise
        error(['requested action unknown: ' opt{iO}.action]);
    end

    if doprint
        print('-deps2', sprintf('plot-%d-%s-%s', iO, analist.adaptid, analist.selectid));
    end
end


function auxPlotCmp(vars, opt, tstr, ind)
if nargin < 4
    C= [vars{opt.nums(1)}.val, vars{opt.nums(2)}.val]';
else
    C= [vars{opt.nums(1)}.val(ind), vars{opt.nums(2)}.val(ind)]';
end
global J
J= C;
%C
nC= size(C,2);
if nC> 0
    plot(C(1,:), C(2,:), 'ok', 'MarkerFaceColor', [0 0 0]);
    hold on
    plot([min(C(1,:)) max(C(1,:))], [0 0], 'k-');
    plot([0 0],[min(C(2,:)) max(C(2,:))],  'k-');
    %    plot(-C(1,:), -C(2,:), 'ok');
    %    plot([min([C(1,:) -C(1,:)]) max([C(1,:) -C(1,:)])], [0 0], 'k-');
    %    plot([0 0],[min([C(2,:) -C(2,:)]) max([C(2,:) -C(2,:)])],  'k-');
    axis tight
    %axis square

    [B,BINT,R,RINT,STATS]= regress(C(2,:)',[C(1,:)' ones(nC,1)]);
    tstr= sprintf('%s: y= %.4f*x %+.4f, \n R^2= %.2f, p= %.4f\n', ...
        tstr, B(1), B(2), STATS(1), STATS(3));
    title(tstr);
    fprintf(1, '%s', tstr);
    vars{opt.nums(1)}.name
    vars{opt.nums(2)}.name
    xlabel(vars{opt.nums(1)}.name);
    ylabel(vars{opt.nums(2)}.name);
end

function out= auxHist(vars, opt, tstr, ind)
if nargin < 4
    C= [vars{opt.nums(1)}.val, vars{opt.nums(2)}.val]';
else
    C= [vars{opt.nums(1)}.val(ind), vars{opt.nums(2)}.val(ind)]';
end
nC= size(C,2);
minc= min(min(C)); maxc= max(max(C));
edge= max([abs(minc), abs(maxc)])*1.01;
borders= linspace(-edge, edge, 11); % n should be odd to separate neg from pos
h= histc(C', borders);
pb= borders(1:end-1)+(borders(2)-borders(1))/2;
bar(pb, h(1:end-1, :));
title(tstr)

m= mean(C');
sd= std(C');
se= sd/sqrt(nC);
for i=1:length(m);
    fprintf(1, '%s: %f +/- %f, (%f)\n', tstr, m(i), se(i), sd(i));
end
fprintf(1, '\ttwo-tailed t-test for mean: ');
for i=1:length(m);
    [H, p]= ttest(C(i,:), 0);
    fprintf(1, 'p%d= %g, \t', i, p);
end
fprintf(1, '\n');
if length(m) == 2
    pval= 1-fcdf((sd(2)/ sd(1))^2, nC-1, nC-1);
    fprintf(1, '\tone-tailed F-test, p= %g\n', pval);
    fprintf(1, '\ttwo-tailed F-test, p= %g\n', testVar(C(1,:), C(2,:)));
end

