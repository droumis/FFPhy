% POPLGROW: Models a proportional population increase, predicting the number of
%           years required to read a given threshold.
%
%     Usage:  [years,N] = poplgrow(initial,proportion,threshold)
%
%           initial = initial population size.
%           proportion = proportional yearly increase in size 0-1).
%           threshold = threshold (target) size.
%           -------------------------------------------------------
%           years = number of years to attain the threshold.
%           N = row vector of population sizes.
%

function [years,N] = poplgrow(initial,proportion,threshold)

  if (proportion < 0)
    error('Proportion must be non-negative');
  end;


  years = 0;
  n = initial;
  N = [];
  
  while (n < threshold)
    years = years + 1;
    n = ceil(n + n*proportion);
    N = [N n];
  end;
  
  close all;
  plot(N);

  return;
  