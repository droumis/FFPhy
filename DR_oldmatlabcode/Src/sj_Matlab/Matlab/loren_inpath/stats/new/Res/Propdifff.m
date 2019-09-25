% PROPDIFFF:  Objective function for propdiff(). Returns absolute difference in 
%             proportions between two groups based on binary observations.
%             Assumes group labels are 1,2.
%
%     Usage: d = propdifff(x,g)
%

% RE Strauss, 1/17/00

function d = propdifff(x,g,not_used1,not_used2)
  x1 = x(g==1);
  x2 = x(g==2);
  p1 = sum(x1>0)/length(x1);
  p2 = sum(x2>0)/length(x2);
  d = abs(p1-p2);

  return;
