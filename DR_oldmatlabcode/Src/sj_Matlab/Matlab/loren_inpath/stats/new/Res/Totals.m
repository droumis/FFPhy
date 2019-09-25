% TOTALS: For an input matrix of any size, returns the column totals, row totals,
%         and grand total.
%
%   Usage: [ColumnTotals,RowTotals,GrandTotal] = totals(X)
%
%         X =            [r x c] input matix.
%         -------------------------------------------------
%         ColumnTotals = [1 x c] vector of column totals.
%         RowTotals =    [r x 1] vector of row totals.
%         GrandTotal =   [1 x 1] grand total.
%

% RE Strauss, 9/5/02

function [ColumnTotals,RowTotals,GrandTotal] = totals(X)

  ColumnTotals = sum(X);              % Calculate column totals
  RowTotals = sum(X')';
  GrandTotal = sum(ColumnTotals);

  return;
  


