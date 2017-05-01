
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code takes the output from runfiles.m
%and computes the within curve probability that
%one trial is different from another
%Anne Smith July 2004

%allps is a matrix that contains the final p-value comparisons
%between trial k and trial j


allps = [];

numbermcs = 10000;   %number of Monte Carlo samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%big loop through all the comparisons

for trial1 = 1:length(qnew)

 pvalue = NaN*zeros(1, length(qnew));

 for trial2 = trial1 + 1: length(qnew)

  shortvec = trial1: trial2-1;

  a1       = eye(size(a,1));

% Appendix D
  for ppp = trial2 - 1: -1: trial1
    a1  = mtimes(a(ppp), a1);
  end

  a1   = mtimes(a1, signewsq(trial2));
  covx = a1(1,1);
 
  mean1 = [ qnew(trial1) qnew(trial2)];
  cov1  = [ signewsq(trial1) covx;  ...
            covx signewsq(trial2) ];
 
  r         = mvnrnd(mean1, cov1, numbermcs);
  r1        = r(:,1); 
  r2        = r(:,2);
  
  pp1       = exp(r1);
  pp2       = exp(r2);
 
  pvalue(trial2)  = length(find(pp1>pp2))/ numbermcs;

 end

allps = [allps; pvalue];

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allps2 = allps;

for ii = 1:length(qnew)
 for jj = ii+1:length(qnew)
  allps2(jj,ii) = NaN; %  make the matrix lower triangular
 end
end

figure(1); pno = 2;
subplot(2, 1, pno)
imagesc(allps2,[0 1]); colormap('bone')

hold on;
co = 1;
[i,j] = find(allps2<0.05/co);
%plot(i,j,'.r');
plot(j,i,'.r');
minj = min(j);
if(~isempty(minj))
 title(['Earliest trial signif above estimated start distribution '  num2str(minj(1)) ])
else
  title('No trials above start distribution')  
end
[i,j] = find(allps2<0.025/co);
%plot(i,j,'*r');
plot(j,i,'*r');

[i,j] = find(allps2>0.95);
%plot(i,j,'.b');
plot(j,i,'.b');
[i,j] = find(allps2>0.975);
%plot(i,j,'*b');
plot(j,i,'*b');
axis([0.5 length(qnew)-.5 0.5 length(qnew)-.5]);
line([0.5 length(qnew)-.5],[0.5 length(qnew)-.5]);

xlabel('Trial Number')
ylabel('Trial Number')
