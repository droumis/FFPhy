% MORTGAGE: Calculates the monthly payment and real cost of a mortgage with 
%           simple yearly interest.
%
%     Usage: [per_month,tot_int,tot_cost] = mortgage(loan,rate,years,{noprint})
%
%           loan =      loan amount, in dollars, after down payment.
%           rate =      annual interest rate, as percentage.
%           years =     duration of loan, in years.
%           noprint =   optional boolean flag indicating that printed results
%                         (to screen) are to be suppressed.
%           ---------------------------------------------------------
%           per_month = monthly payment.
%           tot_int =   total interest.
%           tot_cost =  total cost of loan.
%

% RE Strauss, 3/15/98

function [per_month,tot_int,tot_cost] = mortgage(loan,rate,years,noprint)
  if (nargin < 4)
    noprint = [];
  end;
  if (isempty(noprint))
    noprint = 0;
  end;

  if (rate < 1)                   % Convert proportional rate to percentage
    rate = rate * 100;
  end;

  if (~noprint)
    s = sprintf('\nLoan amount =    %10.2f',loan);
    disp(s);
    s = sprintf('Interest rate =     %7.3f',rate);
    disp(s);
    s = sprintf('Loan duration =     %7.0f',years);
    disp(s);
  end;

  rate = rate ./ 12 ./ 100;       % Convert interest to monthly, decimal
  payments = years .* 12;         % Total monthly payments

  p = rate./(1-(1-rate).^payments);
  per_month = loan * p;

  tot_cost = per_month * payments;
  tot_int = tot_cost - loan;

  if (~noprint)
    s = sprintf('Monthly payment =   %7.2f',per_month);
    disp(s);
    s = sprintf('Total interest = %10.2f',tot_int);
    disp(s);
    s = sprintf('Total cost =     %10.2f \n',tot_cost);
    disp(s);
  end;

  return;
