% RENKONEN: Renkonen's measure of dissimilarity between two communities (sites), 
%           based on relative abundances (proportions) of individuals of 
%           composite species.  See Krebs (1989:304-305).
%
%     Syntax: [dist,abund] = renkenon(abund)
%
%           abund = [N x S] matrix of abundance (counts or proportions of 
%                     individuals) for N taxa across S localities.
%           ---------------------------------------------------------------
%           dist =  [S x S] symmetric distance matrix.
%           abund = input abundance matrix converted to proportions.
%

% Renkonen, O. 1938. Statisch-okologische Untersuchungen uber die terrestiche 
%   Kaferwelt der finnischen Bruchmoore.  Ann. Zool. Soc. Bot. Fenn. Vanamo 
%   6:1-231.
% Krebs, C.J. 1989.  Ecological Methodology.  Harper & Row.

% RE Strauss, 2/22/97 (modified from percdist.m)

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
