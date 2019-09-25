% FIOF:   Improvement-of-fit F test (=general linear test).
%           Input values may be scalars or column vectors or some combination of 
%         the two.  If at least one column vector is present, scalars are 
%         expanded to vectors and statistics are returned for each row of input.
%
%     Usage: [F,pr,df] = fiof(sser,ssef,dfr,dff)
%
%         sser = sse for reduced (restricted) model, that having the fewer 
%                  parameters.
%         ssef = sse for full model, that having the greater number of 
%                  parameters (ssef <= sser).
%         dfr =  sse degrees of freedom for reduced model.
%         dff =  sse degrees of freedom for full model (dfr > dff).
%         ----------------------------------------------------------------------
%         F =    test statistic.
%         pr =   probability under the null hypothesis that the restricted model 
%                  holds.
%         df =   2-element row vector containing the degrees of freedom for 
%                  the test (dfr-dff, dff).
%

% Neter J, W Wasserman, MH Kutner. 1985. Applied linear statistical models: 
%   regression, analysis of variance, and experimental designs.  2nd ed.  
%   Richard D. Irwin, Homewood IL.

% RE Strauss, 1/3/01

function [F,pr,df] = fiof(sser,ssef,dfr,dff)
  sser = sser(:);                     % Convert input to col vectors
  ssef = ssef(:);
  dfr = dfr(:);
  dff = dff(:);

  [b,sser,ssef,dfr,dff] = samelength(sser,ssef,dfr,dff);  % Check for same length
  if (~b)
    error('  FIOF: input matrices not consistent in length.');
  end;

  d = dfr-dff;
  i = find(d<eps | ssef<eps | dff<eps);
  if (~isempty(i))
    nf = NaN*ones(length(i),1);
    d(i) = nf;
    ssef(i) = nf;
    dff(i) = nf;
  end;

  F = (sser-ssef)./(dfr-dff) ./ (ssef./dff);
  F = max([F,zeros(length(F),1)]')';  % Convert neg values to zeros
  df = [dfr-dff dff];
  pr = 1 - fcdf(F,df(:,1),df(:,2));

  return;
