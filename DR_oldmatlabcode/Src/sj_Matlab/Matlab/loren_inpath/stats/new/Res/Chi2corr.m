% CHI2CORR: Chi-square 'correlations' among the columns of a contingency table
%           of counts.
%
%     Syntax: corr = chi2corr(table)
%
%           table = [r x c] table of counts
%           dist =  [c x c] symmetric matrix of column correlations
%

% Davis, JC. 1986. Statistics and Data Analysis in Geology, 2nd ed, p. 581. Wiley.

% RE Strauss, 11/27/95

function corr = chi2corr(table)
  [r,c] = size(table);                % Table size

  total = sum(sum(table));            % Grand total of counts
  table = table / total;              % Convert elements to joint probabilities
  coltot = sum(table);                % Marginal totals of probabilities
  rowtot = sum(table');

  for i = 1:r                         % Sqrts of 'chi-square' cell contributions
    for j = 1:c
      exp = rowtot(i)*coltot(j);
      table(i,j) = (table(i,j)-exp)/sqrt(exp);
    end;
  end;

  corr = table' * table;              % Correlations among columns

  return;
