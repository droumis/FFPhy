function auxCmpVars(var1, var2, opt)
%function auxCmpVars(vars1, vars2, opt)
% 
%  Show relationship between two variables

if nargin < 3; opt= []; end

if isfield(opt, 'varid') & opt.varid
    v1= var1.val(:,opt.varid); 
    v2= var2.val(:,opt.varid); 
    n1= length(v1);
    n2= length(v2);
else
    n1= prod(size(var1.val)); 
    n2= prod(size(var2.val)); 
    v1= reshape(var1.val,n1,1); 
    v2= reshape(var2.val,n2,1);
end

%keyboard

showtitle=0; 
if isfield(var1, 'name') showtitle= 1; end
bothtitles=0; 
if isfield(var1, 'name') & isfield(var2, 'name') & strcmp(var1.name, var2.name)
    bothtitles=1;
end

C= [v1, v2]';
%C
nC= size(C,2);
if nC> 0
    plot(C(1,:), C(2,:), 'ok', 'MarkerSize', 4, 'MarkerFaceColor', [0 0 0]);
%    keyboard
    hold on
    plot([min(C(1,:)) max(C(1,:))], [0 0], 'k-');
    plot([0 0],[min(C(2,:)) max(C(2,:))],  'k-');
    %    plot(-C(1,:), -C(2,:), 'ok');
    %    plot([min([C(1,:) -C(1,:)]) max([C(1,:) -C(1,:)])], [0 0], 'k-');
    %    plot([0 0],[min([C(2,:) -C(2,:)]) max([C(2,:) -C(2,:)])],  'k-');
    axis tight
    %axis square

    [B,BINT,R,RINT,STATS]= regress(C(2,:)',[C(1,:)' ones(nC,1)]);
    if showtitle & ~bothtitles
        tstr= sprintf('%s: y= %.4f*x %+.4f, \n R^2= %.2f, p= %.4f\n', ...
            var1.title, B(1), B(2), STATS(1), STATS(3));
    elseif  bothtitles
        tstr= sprintf('%s, %s: y= %.4f*x %+.4f, \n R^2= %.2f, p= %.4f\n', ...
            var1.title, var2.title, B(1), B(2), STATS(1), STATS(3));
    end
    title(tstr);
    fprintf(1, '%s', tstr);
    xlabel(var1.name);
    ylabel(var2.name);
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

