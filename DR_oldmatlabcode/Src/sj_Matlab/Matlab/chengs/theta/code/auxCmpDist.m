function out= auxCmpDist(var1, var2, opt)
%function out= auxCmpDist(var1, var2, opt)
% 
%  Compare two samples of the same random variable
% 

if nargin < 3; opt= []; end

if ~isfield(opt, 'varid') | ~opt.varid
    opt.varid= 1;
end
v1= var1.val(:,opt.varid); 
v2= var2.val(:,opt.varid); 
n1= length(v1);
n2= length(v2);

%if opt.varid== 2; keyboard; end

bothname= 0;
if ~strcmp(var1.name, var2.name); 
    warning('comparing two different variables'); 
    bothname= 1;
end

plottype='pdf'; if isfield(opt, 'plot') plottype= opt.plot; end
lw=3; if isfield(opt, 'LineWidth') lw= opt.LineWidth; end

if strcmp(plottype, 'cdf')
    % plot cdf
    s1= sort(v1); y1= [1:n1]'/n1;
    s2= sort(v2); y2= [1:n2]'/n2;
    plot(s1, y1, 'LineWidth', lw); hold on
    plot(s2, y2, 'r', 'LineWidth', lw); hold on
    ylabel('cdf')
%    keyboard
else

    nbins= 10; if isfield(opt, 'nbins') nbins= opt.nbins; end


    minv= min([min(v1), min(v2)]); 
    maxv= max([max(v1), max(v2)]); 
    %edge= max([abs(minv), abs(maxv)])*1.01;
    %borders= linspace(-edge, edge, 11); % n should be odd to separate neg from pos
    borders= linspace(.95*minv, 1.05*maxv, nbins+1);
    pb= borders(1:end-1)+(borders(2)-borders(1))/2;
    h1= histc(v1, borders)/ n1;
    h2= histc(v2, borders)/ n2;
    if strcmp(plottype, 'pdf')
        plot(pb, h1(1:end-1), 'LineWidth', lw); hold on
        plot(pb, h2(1:end-1), 'r', 'LineWidth', lw);
        ylabel('pdf')
    elseif strcmp(plottype, 'hist')
%        h1= histc(v1, borders);
%        h2= histc(v2, borders);
        bar(pb, [h1(1:end-1), h2(1:end-1)]);
        ylabel('count')
    end
end

if ~isfield(opt,'varname') opt.varname=''; end
xlabel([opt.varname var1.name]);
if ~strcmp(var1.name, var2.name)
    warning('compared variables have different names');
end

x= [min([min(v1),min(v2)])-1; sort(v1); max([max(v1),max(v2)])+1;];
y= [0;[1:n1]'/n1; 1];

%[h,p]= kstest(v2, [x, y]);
[h,pKS]= kstest2(v1,v2);
pW= ranksum(v1,v2);

%disp([opt.varname var1.name]);
%if h
%    fprintf(1, 'KS test: %s diff''t from %s, p= %.2g\n', var1.title, var2.title, pKS);
%else
%    fprintf(1, 'KS test: %s NOT diff''t from %s, p= %.2g\n', var1.title, var2.title, pKS);
%end
%fprintf(1, '\tranksum test: p= %.2g\n', pW);


if isfield(opt, 'legend') & opt.legend
    legend({var1.title, var2.title},0);
end

titlestr= '';
if isfield(opt, 'title') titlestr= opt.title; end
if ~isfield(opt, 'testsummary') | opt.testsummary
    if ~isempty(titlestr) titlestr= [titlestr '\n']; end
    if bothname
        titlestr= sprintf('%s%s, %s: KS= %.2g, W= %.2g', titlestr, var1.name, var2.name, pKS, pW);
    else
        titlestr= sprintf('%s%s: KS= %.2g, W= %.2g', titlestr, var1.name, pKS, pW);
    end
end
%        keyboard
title(titlestr);
th=get(gca, 'Title');
set(th, 'Interpreter', 'none')


%m= mean(C');
%sd= std(C');
%se= sd/sqrt(nC);
%for i=1:length(m);
%    fprintf(1, '%s: %f +/- %f, (%f)\n', tstr, m(i), se(i), sd(i));
%end
%fprintf(1, '\ttwo-tailed t-test for mean: ');
%for i=1:length(m);
%    [H, p]= ttest(C(i,:), 0);
%    fprintf(1, 'p%d= %g, \t', i, p);
%end
%fprintf(1, '\n');
%if length(m) == 2
%    pval= 1-fcdf((sd(2)/ sd(1))^2, nC-1, nC-1);
%    fprintf(1, '\tone-tailed F-test, p= %g\n', pval);
%    fprintf(1, '\ttwo-tailed F-test, p= %g\n', testVar(C(1,:), C(2,:)));
%end

%keyboard
