function simCoincDist(N,n1,n2)
%function simCoincDist(N,n1,n2)

%N= 10;
%n1=1;
%n2=1;

reps= 1e3;

x1= zeros(N,1); x2= x1;
x1(1:n1)=1; x2(1:n2)=1;
n12= zeros(reps, 1);
p= zeros(reps, 1);

den= nchoosek(N,n1);

for ir=1:reps
    x1= perm(x1);
    x2= perm(x2);

    n12(ir)= sum(x1&x2);
%    for n=n12(ir):min(n1,n2)
%        p(ir)= p(ir)+nchoosek(n2,n)*nchoosek(N-n2,n1-n)/ den;
%    end
    p(ir)=  nchoosek(n2, n12(ir))*nchoosek(N-n2,n1-n12(ir))/ den;
end

p12= zeros(n1+1, 1);
for n=max(0,n1+n2-N):min(n1,n2)
    p12(n+1)= nchoosek(n2,n)*nchoosek(N-n2,n1-n)/ den;
end

n=[0:n1]';
h12=histc(n12,n);
h12= h12/sum(h12);
figure(1); clf
bar(n, h12);
hold on
plot(n,p12, 'r.')

mu= mean(n12);
V= var(n12);
fprintf(1, '%5s: %5s %5s\n', 'var', 'empir', 'theor');
fprintf(1, '%5s: %5.1f %5.1f\n', 'mean', mu, n1*n2/N);
fprintf(1, '%5s: %5.2f %5.2f\n', 'var', V,  n1*n2*(N-n1)*(N-n2)/(N^2*(N-1)));
np= normpdf(n, mu, sqrt(V));
plot(n, np, 'b--');
np= normpdf(n, mu-0.5, sqrt(V));
plot(n, np, 'g--');

%figure(2); clf
%cdfplot(p)
%maxp= max(p)+1e-15;
%x= linspace(0, maxp, 20)';
%y= (x/ maxp).^2;
%cdfplot(p)
%hold on
%plot(x,y, 'r')
%[rej, pval]= kstest(p, [x y])


