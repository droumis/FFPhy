% CORRPERM: Find the optimally permuted matrix by enumeration.
%
%     Usage: [permute,obj,rc,stats,distrib] = corrperm(r,{nbest})
%
%           r =       square symmetric matrix
%           nbest =   if positive: indicates the number of best permutations
%                       to be saved (default=1)
%                     if negative: indicates a restart from a previously
%                       incomplete run
%           permute = best permutation of matrix
%           obj =     value of objective function
%           rc =      permuted matrix
%           stats =   [improvements x 2] matrix of accumulated performance
%                       statistics:
%                         cumulative total permutations
%                         best accumulative objective-function value
%           distrib = sampling distribution of objective-function values
%
%     Creates three files: 'best.out', containing a record of the n-best
%     permutations; 'iterlog.out', containing a record of the iterations during
%     which improvements in the objective function occur or at breakpoints; and
%     'distrib.mat', a binary file containing a binned distribution of
%     objective-fn values.  These files can be used to restart the procedure
%     (picking up at the most recent improvement or breakpoint) if an interruption
%     occurs.
%

% RE Strauss, 3/18/95

function [permute,obj,rc,stats,distrib] = corrperm(r,nbest)
  cd '\matlab\output';                % Set current directory

  TRUE=1; FALSE=0;
  list_update = FALSE;                % Update files during report
  restart = FALSE;                    % Restart from previous run
  get_stats = FALSE;
  get_distrib = FALSE;

  default_nbest = 1;
  to_report = 10000;                  % Number of iterations till next report

  distrib_min = -0.1;                 % Min interval in distribution
  distrib_max =  0.5;                 % Max interval in distribution
  distrib_delta = 0.001;              % Interval width
  distrib_length = ((distrib_max-distrib_min)/distrib_delta)+1;

  if (nargout > 3)
    get_stats = TRUE;
  end;
  if (nargout > 4)
    get_distrib = TRUE;
  end;

  report_counter = 0;
  iter = 0;                           % Iteration number

  N = size(r,1);
  nfact = prod(1:N);                  % Number of permutations
  w = corradj(N);                     % Adjacency matrix

  if (nargin<2)
    nbest = default_nbest;
  else
    if (nbest < 0)                    % If negative, restart from previous run
      restart = TRUE;                 %   flag restart

      load distrib;                   %   Get distribution
      distrib_length = size(distrib,1);

      load best.out;                  %   Get 'best' matrix
      nbest = size(nbest,1);          %   Determine 'nbest' from number of rows

      load iterlog.out;               %   Get most recent results
      [rp,cp] = size(iterlog);
      iter = iterlog(rp,1);           %   Last iteration
      prev_best = iterlog(rp,2);      %   Last best obj-fn value
      p = iterlog(rp,3:cp);           %   Last permutation
    end;
  end;

  if (~restart)                       % If new analysis, initialize
    p = [N:-1:1];                     %   last permutation (wraps to first)
    best = zeros(nbest,N+2);          %   list of best solutions
    prev_best = 0;                    %   best obj-fn value
    distrib = zeros(distrib_length,1);%   distribution of obj-fn values

    if (exist('e:\matlab\output\iterlog.out')) % Delete previous report log
      delete e:\matlab\output\iterlog.out;
    end;
  end;

  % Iterate thru permutations

  while (iter < nfact)
    iter = iter+1;
    p = permnext(p);                  % Get permutation
    rc = r(p,p);                      % Permute matrix
    obj = corrobj(rc,w);              % Objective function for permuted matrix

    % Index into distribution bin
    d = round((distrib_length-1)*(obj-distrib_min)/(distrib_max-distrib_min))+1;
    distrib(d) = distrib(d)+1;        % Increment bin

    if (obj > best(1,1))              % If best so far, display
      disp(sprintf('  Permutations=%7.0f, best=%7.4f', iter, obj));
      corrprm2(iter,best(1,1),p,distrib);   % Write to files
      best = addlist(best,[obj iter p],0);  % Add to list of nbest
      corrprm1(best);
    else
      % Add to list of nbest if good enough
      [best,update] = addlist(best,[obj iter p],0);
      if (update)
        list_update = TRUE;
      end;
    end;

    report_counter = report_counter+1;
    if (report_counter == to_report)  % If at milestone, display
      disp(sprintf('  Permutations=%7.0f, best=%7.4f', iter, best(1,1)));
      corrprm2(iter,best(1,1),p,distrib);

      if (list_update)
        corrprm1(best);               % Update current results to file
        list_update = FALSE;
      end;
      report_counter = 0;             % Reset counter
    end;
  end; % permutations

  obj = best(:,1);                    % Vector of best obj-fn values
  permute = best(:,3:N+2);            % Matrix of best permutations
  p = permute(1,:);                   % Best permutation
  rc = r(p,p);                        % Best matrix

  for i=2:nbest                       % Chk congruency of best permutations
    forward = 0;
    for j=1:N                         % Chk pos of each element against
      k = find(p == permute(i,j));    %   pos in best permutation
      forward = forward + (j-k)^2;    % Accum squared differences in position
    end;
    backward = 0;
    for jj=N:-1:1                     % Repeat for reverse permutation
      j = N-jj+1;
      k = find(p == permute(i,jj));
      backward = backward + (j-k)^2;
    end;
    if (backward < forward)           % Reverse permutation if necessary
      permute(i,1:N) = permute(i,N:-1:1);
    end;
  end;

  disp(sprintf('  Permutations=%7.0f, best=%7.4f', iter, best(1,1)));
  corrprm1(best);                     % Write final results to files
  corrprm2(iter,best(1,1),p,distrib);

  % Recover stats from file if requested
  if (get_stats)
    load iterlog.out;
    ri = size(iterlog,1);
    stats = iterlog(1,[1:2]);
    best = stats(1,2);

    for i=2:(ri-1)
      if (iterlog(i,2) ~= best)
        iter = iterlog(i,1);
        stats = [stats; iter-1, best];
        best = iterlog(i,2);
        stats = [stats; iter, best];
      end;
    end;
    stats = [stats; iterlog(ri,[1:2])];

    figure(1);
    plot(stats(:,1),stats(:,2));
  end;

  % Recover distribution from file if requested
  if (get_distrib)
    load distrib;
    d = [distrib_min : distrib_delta : distrib_max]';
    indx = find(distrib);
    first = indx(1);
    last = indx(length(indx));
    distrib = [d(first:last) distrib(first:last)/2];

    figure(2);
    bar(distrib(:,1),distrib(:,2));
  end;

  return;
