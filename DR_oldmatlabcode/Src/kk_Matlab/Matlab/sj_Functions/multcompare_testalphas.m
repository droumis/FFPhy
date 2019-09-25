function [cmp,whicha] = multcompare_testalphas(stats, dimension, varargin)

possible_alphas = [0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
[otherArgs] = procOptions(varargin);

[c] = multcompare(stats,'dimension',dimension,'display','off');
N = size(c,1);

cmp = nan(N,3);
whicha = zeros(N,1);

for a = 1:length(possible_alphas)
  cmp_{a} = multcompare(stats,'dimension',dimension,'alpha',possible_alphas(a),'display','off');
  flag = 0;
  for i = 1:N
    cmp(i,1:2) = cmp_{a}(i,1:2);
    if ((cmp_{a}(i,3) > 0) == (cmp_{a}(i,5) > 0))
      flag = 1;
      cmp(i,3) = possible_alphas(a);
      whicha(i) = a;
    end
  end
  if (flag == 0)
    break;
  end
end



