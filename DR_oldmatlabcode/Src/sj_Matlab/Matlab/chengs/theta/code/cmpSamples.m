function cmpSamples(x,y,alpha)
%function cmpSamples(x,y,alpha)

x= x(find(isfinite(x)));
y= y(find(isfinite(y)));

if nargin < 3; alpha= 0.05; end

fprintf(1, 'Two-sample Kolmogorov-Smirnov\n');
[h,p]= kstest2(x,y,alpha,0);
fprintf(1, '\tH0: F(x)=F(y), H1: F(x) ~= F(y), p= %.4f\n', p);
[h,p]= kstest2(x,y,alpha,1);
fprintf(1, '\tH0: F(x)=F(y), H1: F(x)  > F(y), p= %.4f\n', p);
[h,p]= kstest2(x,y,alpha,-1);
fprintf(1, '\tH0: F(x)=F(y), H1: F(x)  < F(y), p= %.4f\n', p);

fprintf(1, 'Wilcoxon rank sum\n');
[p,h]= ranksum(x,y,alpha);
fprintf(1, '\tH0: median(x)=median(y), p= %.4f\n', p);


[h, p]= ttest2(x,y, alpha, 0);
fprintf(1, 'two sample t-test F(x)=F(y), p=%.4g\n', p);
