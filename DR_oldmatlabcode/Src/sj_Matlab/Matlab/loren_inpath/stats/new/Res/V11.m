% V11:  Estimate V(n), the variance of the largest order statistic, using the method of
%       Shea & Scallon (1988; AS 128).  This is needed for the Davis & Stephens (1978)
%       algorithm for approximating the variance-covariance matrix of ordered sample of
%       independent observations from a standard normal distribution.
%
%     Usage: v = v11(n)
%

% Shea, BL and AJ Scallon. 1988. Remark AS R72: a remark on algorithm AS 128. Approximating
%   the covariance matrix of normal order statistics.  Appl. Stat. 38:151-155.

% RE Strauss, 12/9/02

function v = v11(n)

  mpt15 = -0.15;
  pt09 = 0.091105452691946d0;
  a0 = 0.04619831847696d0;
  a1 = -0.147930264017706d0;
  a2 = -0.451288155800301d0;
  a3 = 0.010055707621709d0;
  a4 = 0.007412441980877d0;
  a5 = -0.001143407259055d0;
  a6 = 0.54428754576d-4;
  b0 = -0.934d-4;
	b1 = -0.5950321d0;
  b2 = 0.0165504d0;
  b3 = 0.0056975d0;
  b4 = -0.8531d-3;
  c0 = 0.7956d-11;
  c1 = -0.595628869836878d0;
  c2 = 0.08967827948053d0;
  c3 = -0.007850066416039d0;
  c4 = -0.296537314353d-3;
  c5 = 0.215480033104d-3;
  c6 = -0.33811291323d-4;
  c7 = 0.2738431187d-5;
  c8 = -0.106432868d-6;
  c9 = 0.1100251d-8;
  d0 = 0.093256818332708d0;
  d1 = 1.336952989217635d0;
  d2 = -1.783195691545387d0;
  d3 = 0.488682076188729d0;
  d4 = -0.078737246197474d0;
  d5 = 0.00662561987806d0;
  d6 = -0.226486218258d-3;

	v = 1;
  if (n<1)
    error('  V11: n < 1.');
  elseif (n==1)
    return;
  end;
  
	x = n;
	if (n>370) 
	  x = (x^mpt15-1)/mpt15;
	  v = exp(b0+x*(b1+x*(b2+x*(b3+x*b4))));
	elseif (n<=100) 
	  x = (x^pt09-1)/pt09;
	  v = exp(c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8+x*c9)))))))));
	elseif (n<=200) 
	  x = log(a0 + x);
	  v = exp(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + x*a6)))));
	else
	  x = log(d0 + x);
	  v = exp(d1 + x*(d2 + x*(d3 + x*(d4 + x*(d5 + x*d6)))));
	end;

	return;
