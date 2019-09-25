% STEPDISC: Stepwise discriminant analysis to determine the best subset of variables.  
%           Introduces variables so as to maximize at each step the Lawley-Hotelling 
%           trace (=Rao's V).  This measure is proportional to the mean Mahalanobis 
%           distance.  Optionally boostraps the data matrix and returns an accumulated 
%           matrix of sequences of introduced variables as well as matrices of diagonstics.
%             If both group and subgroup identifiers are provided for each observation 
%           (e.g., identifying species and population within species), discriminant 
%           analyses are done independently by group at each step (discriminating the 
%           subgroups) and the resultant values of V averaged.  This procedure allows 
%           the identification of common sets of variables that best discriminate 
%           subgroups.  Such sets are likely to be suboptimal for any particular group, 
%           but are optimal (by the mean-V criterion) across groups.
%
%     Usage: [incl,Vcum,Vpcum,bootperc,bootbest,bootplace,bootincl,bootVpcum] = ...
%                 stepdisc(X,grps,{maxvars},{iter})
%
%         X =        [n x p] data matrix (obs x vars).
%         grps =     [n x 1] vector of group identifiers, or [n x 2] matrix of group 
%                      (col 1) and subgroup (col 2) identifiers. 
%         maxvars =  maximum number of variables to be included in model [default = all 
%                      variables, or the maximum supported by the observations].
%         iter =     optional number of bootstrap iterations [default = 0].
%         ------------------------------------------------------------------------------
%         incl =     list of indices of variables in order of inclusion.
%         Vcum =     maximum V at each step.
%         Vpcum =    maximum Vp (=V divided by number of variables) at each step.
%         bootperc = percentage of bootstrap iterations in which each variable in 'incl' 
%                      is in the corresponding position (step).
%         bootbest = [step x nvars] matrix giving the best subset of variables for each
%                      step, based on bootstrap results.
%         bootplace = [p x maxvars] matrix giving percentage of bootstrap iterations in 
%                      which variable p_i is in step maxvars_i.
%         bootincl = matrix of lists of variables, one row per bootstrap iteration.
%         bootVcum = vector of maxima of the Vpcum function.
%

% RE Strauss, 7/30/98
%   11/29/99 - changed calling sequence.
%    6/16/00 - added check for missing data.

function [incl,Vcum,Vpcum,bootperc,bootbest,bootplace,bootincl,bootVpcum] = ...
            stepdisc(X,grps,maxvars,iter)
  
  if (nargin < 3) maxvars = []; end;
  if (nargin < 4) iter = []; end;

  [n,p] = size(X);
  [r,c] = size(grps);

  supergrps = [];
  if (r==1 & c>1)                   % If group vector is a row, transpose to column
    grps = grps';
  elseif (r>1 & c==2)               % If two columns are passed,
    supergrps = grps(:,1);            % col 1 = grps
    grps = grps(:,2);                 % col 2 = subgrps
  end;

  if (misscheck(X,grps))
    error('  STEPDISC: data matrix or grouping vector contains missing data.');
  end;

  supergrp_ids = uniquef(supergrps);
  nsgrps = length(supergrp_ids);

  if (isempty(maxvars))             % Default argument values
    maxvars = p;
  end;
  if (isempty(iter))
    iter = 0;
    bootperc = [];
    bootplace = [];
    bootincl = [];
  end;
                                    % Allocate output matrices
  vars_incl = zeros(1,p);             % Positions in which vars introduced
  Vcum =  zeros(1,maxvars);           % Cumulative V
  Vpcum = zeros(1,maxvars);           % Cumulative V per variable

  if (~isempty(supergrps))
    Vg =  zeros(nsgrps,1);
    Vpg = zeros(nsgrps,1);
  end;

  for step = 1:maxvars              % Add vars one at a time
    vi = find(vars_incl>0);           % Indices of vars included
    vni = find(vars_incl==0);         % Indices of vars not included
    nvni = length(vni);               % Number of vars not included

    V =  zeros(nvni,1);               % Get objective fn value for each unused var
    Vp = zeros(nvni,1);
    for v = 1:nvni                    % Cycle thru unused vars
      if (isempty(supergrps))           % Get values for V & Vp
        [V(v),Vp(v)] = lawley(X(:,[vi vni(v)]),grps);
      else
        for g = 1:nsgrps                % If supergrps present, cycle thru them
          i = find(supergrps==g);         % Isolate observations in current supergrp
          Xg = X(i,:);                    % Get data
          [Vg(g),Vpg(g)] = lawley(Xg(:,[vi vni(v)]),grps(i));
        end;
        V(v) =  mean(Vg);                 % Average stats across supergrps
        Vp(v) = mean(Vpg);
      end;
    end;

    if (finite(sum(V)))               % If all fn values are finite,
      [Vpmax,i] = max(Vp);            %   select largest
      vars_incl(vni(i(1))) = step;
      Vpcum(step) = Vpmax;
      Vcum(step) = V(i);
    else                              % Else quit here
      break;
    end;
  end;

  [y,incl] = sort(vars_incl);
  i = find(y>0);
  incl = incl(i);
  len_incl = length(incl);
  Vcum = Vcum(1:len_incl);
  Vpcum = Vpcum(1:len_incl);


  if (iter)                             % Bootstrap results
    bootincl =  zeros(iter,maxvars);
    bootperc =  zeros(size(incl));
    bootplace = zeros(p,maxvars);
    bootbest =  zeros(maxvars,maxvars);
    bootVpcum = zeros(iter);

    for it = 1:iter
%it
      Xb = bootsamp(X,grps);
      vars_incl = zeros(1,p);

      for step = 1:maxvars
        vi = find(vars_incl>0);             % Indices of vars included
        vni = find(vars_incl==0);           % Indices of vars not included
        nvni = length(vni);                 % Number of vars not included

        V = zeros(nvni,1);
        Vp = V;
        for v = 1:nvni
          if (isempty(supergrps))           % Get values for V & Vp
            [V(v),Vp(v)] = lawley(Xb(:,[vi vni(v)]),grps);
          else
            for g = 1:nsgrps                % If supergrps present, cycle thru them
              i = find(supergrps==g);         % Isolate observations in current supergrp
              Xbg = X(i,:);                   % Get data
              [Vg(g),Vpg(g)] = lawley(grps(i),Xbg(:,[vi vni(v)]));
            end;
            V(v) =  mean(Vg);                 % Average stats across supergrps
            Vp(v) = mean(Vpg);
          end;
        end;

        if (finite(sum(V)))
          [Vpmax,i] = max(Vp);
          bootVpcum(it) = Vpmax;
          vars_incl(vni(i(1))) = step;
        else
          break;
        end;
      end;

      [y,inc] = sort(vars_incl);
      i = find(y>0);
      bootincl(it,:) = [inc(i) zeros(1,maxvars-length(i))];
    end;

    % Find the percentage of times (bootstrap iterations) that each variable is in it's
    % given step position.

    for i = 1:length(incl)
      bootperc(i) = 100*sum(bootincl(:,i)==incl(i))/iter;
    end;
    for i = 1:p
      for j = 1:maxvars
        bootplace(i,j) = 100*sum(bootincl(:,j)==i)/iter;
      end;
    end;

    % Accumulate the sums (cumsum) of rows in matrix 'bootplace'.  For each given number of
    % variables (step, =p_i), select the p_i variables having the largest cum-sum 
    % representation.

    bootcum = cumsum(bootplace')';
bootcum
    for i = 1:maxvars
      [b,j] = sort(-bootcum(:,i));
      bootbest(i,1:i) = j(1:i)';
    end;
  end;

  return;
