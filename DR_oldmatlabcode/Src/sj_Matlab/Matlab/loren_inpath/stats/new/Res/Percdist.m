% PERCDIST: The distance measure corresponding to the measure of 
%           percentage similarity.  See renkonen().
%
%     Syntax: [dist,abund] = percdist(abund)
%
%           abund = [N x S] matrix of abundance (counts or proportions of 
%                     individuals) for N taxa across S localities.
%           --------------------------------------------------------------
%           dist =  [S x S] symmetric distance matrix.
%           abund = input abundance matrix converted to proportions.
%

% RE Strauss, 2/22/97

function [dist,abund] = percdist(abund)
  [n,s] = size(abund);
  sim = ones(s,s);

  sumabund = sum(abund);            % Convert to proportions by column
  for j = 1:s
    if (sumabund(j)>0)
      abund(:,j) = abund(:,j)/sumabund(j);
    end;  
  end;

  for i = 1:(s-1)                   % All possible pairs of sites (columns)
    for j = (i+1):s
      sim(i,j) = sum(min([abund(:,[i,j])]'));
      sim(j,i) = sim(i,j);
    end;
  end;

  dist = 1-sim;                     % Convert to distance matrix

  return;
