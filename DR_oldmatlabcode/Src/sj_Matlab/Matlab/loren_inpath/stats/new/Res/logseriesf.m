% LOGSERIESF: objective function
%
%     Usage: sse = logseriesf([alpha,x],obs_abund)

function sse = logseriesf(params,obs_abund)
  alpha = params(1);
  x = params(2);

  nspecies = length(obs_abund);
  max_ind = max(obs_abund);
  nind = [1:max_ind]';
  ffind = zeros(size(nind));

  [u,f] = uniquef(obs_abund,1);
  for i = 1:length(u)
    j = find(nind==u(i));
    if (~isempty(j))
      ffind(j) = f(i);
    end;
  end;
  obs_abund = ffind;

  exp_abund = zeros(max_ind,1);
  for i = 1:max_ind
    exp_abund(i) = alpha * (x.^i)/i;
  end;
%  sum_exp = sum(exp_abund);
%  exp_abund = nspecies*exp_abund/sum_exp;

  d = obs_abund - exp_abund;
  sse = d'*d;
 
  return;

  