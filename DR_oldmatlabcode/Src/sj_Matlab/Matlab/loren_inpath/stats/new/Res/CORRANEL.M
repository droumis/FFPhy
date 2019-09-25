% CORRANEL: Find the optimally permuted matrix by simulated annealing.
%
%     Syntax: [permute,obj,rc,stats,distrib] = corranel(r,nbest)
%
%           r =       square symmetric matrix
%           nbest =   saves the n-best permutations (default=1)
%           ------------------------------------------------------------
%           permute = best permutation of matrix
%           obj =     value of objective function
%           rc =      permuted matrix
%           stats =   [steps x 4] matrix of accumulated performance
%                       statistics:
%                         temperature step
%                         objective-function value at end of step
%                         best accumulative objective-function value
%                         cumulative total permutations
%           distrib = sampling distribution of objective-function values
%

% RE Strauss, 12/20/95

function [permute,obj,rc,stats,distrib] = corranel(r,nbest)
  if (nargin < 2) nbest = []; end;

  TRUE = 1; FALSE = 0;
  RANDOM = 1; SWAP = 2; INVERT = 3; SPLICE = 4;
  default_nbest = 1;
  default_distrib = 1000;

  if (isempty(nbest))
    nbest = default_nbest;
  end;

  N = size(r,1);                  % Size of square symmetric matrix
  max_factor = 1.5;
  max_tsteps = 100*max_factor;    % Max number of temperature steps
  max_pertemp = 100*N*max_factor; % Max permutations per temperature step
  max_success = 10*N*max_factor;  % Max successes before next temperature
  init_temp = .5;                 % Initial temperature
  tp = 0.65;                      % Target proportion of temperature steps
                                  %   for 1 replacement per step

  distrib_length = default_distrib;
  distrib = zeros(1,distrib_length);  % Sampling distrib of obj-fn values

  % Annealing schedule reduction factor
  temp_reduce = exp((-(log(init_temp)+log(max_pertemp)))/(tp*max_tsteps));
disp(temp_reduce);

  clk = clock;                    % Seed random-number generator
  rand('seed',clk(6));            %   from system clock

  w = corradj(N);                 % Adjacency matrix
  best = zeros(nbest,N+1);        % List of best solutions
  permute = 1:N;                  % Initial permutation
  obj = corrobj(r,w);             % Initial objective-function value
  best(1,:) = [obj permute];      % Current solution = best so far

  get_stats = FALSE;
  if (nargout > 3)
    get_stats = TRUE;
    stats = [0 obj best(1,1) 0];    % Accumulated statistics
  end;

  temp = init_temp;
  total_perms = 0;

  for t=1:max_tsteps              % Cycle thru temperatures
    nsuccess = 0;
    for p=1:max_pertemp           % Try permutations
      total_perms = total_perms+1;
%      kindperm = ceil(rand*4);    % Determine kind of permutation
      kindperm = 1;               % Random only

      if (kindperm == RANDOM)     % Completely random permutation
        new_permute = randperm(N);
      end;  % Random

      if (kindperm == SWAP)       % Swap two elements
        n1 = ceil(rand*N);        % Choose first
        n2 = ceil(rand*(N-1));    % Choose second
        if (n2 >= n1)
          n2 = n2+1;
        end;
        new_permute = permute;
        xx = new_permute(n1);
        new_permute(n1) = new_permute(n2);
        new_permute(n2) = xx;
      end;  % Swap

      if (kindperm == INVERT)     % Invert the segment
        nn = 0;
        while (nn<3)                          % Find segment
          n1 = ceil(rand*N);                  % Beginning
          n2 = wrap(n1+ceil(rand*(N-4)),N);   % End
          nn = rem((n1-n2+N-1),N);            % Number of elements not in segment
        end;

        i1 = n1-1;                            % Invert sequence, wrapping around
        i2 = n2+1;                            %   if necessary
        new_permute = permute;
        for i=1:floor(N-nn)/2
          i1 = wrap(i1+1,N);
          i2 = wrap(i2-1,N);
          xx = new_permute(i1);
          new_permute(i1) = new_permute(i2);
          new_permute(i2) = xx;
        end;
      end;  % Inversion

      if (kindperm == SPLICE)     % Splice the segment
        nn = 0;
        while (nn<3)                          % Find segment
          n1 = ceil(rand*N);                  % Beginning
          n2 = wrap(n1+ceil(rand*(N-4)),N);   % End
          nn = rem((n1-n2+N-1),N);            % Number of elements not in segment
        end;

        n3 = wrap((n2+ceil((nn-1)*rand)),N);  % Choose target location
        n4 = wrap(n3+1,N);                    % Element following n3
        n5 = wrap(n1-1,N);                    % Element preceding n1
        n6 = wrap(n2+1,N);                    % Element following n2
        m1 = 1+rem((n2-n1+N),N);              % Number from n1-n2
        m2 = 1+rem((n5-n4+N),N);              % Number from n4-n5
        m3 = 1+rem((n3-n6+N),N);              % Number from n6-n3

        temp_permute = zeros(size(permute));
        p=0;
        for i=1:m1                            % Copy chosen segment n1-n2
          p=p+1;
          temp_permute(p) = permute(wrap(n1+i-1,N));
        end;
        for i=1:m2                            % Copy segment n4-n5
          p=p+1;
          temp_permute(p) = permute(wrap(n4+i-1,N));
        end;
        for i=1:m3                            % Copy segment n6-n3
          p=p+1;
          temp_permute(p) = permute(wrap(n6+i-1,N));
        end;

        if (any([n2<n1, n1==1]))   % If position 1 was in segment
          ls = m1+m2;              %   left-shift n6 back to position 1
        else;                      % Else left-shift original position 1
          for i=1:N                %   back to position 1
            if (temp_permute(i)==permute(1))
              ls = i-1;
              break;
            end;
          end;
        end;
        for i=1:N                   % Left-shift
          new_permute(i) = temp_permute(wrap(i+ls,N));
        end;
      end;  % Splice

      % If new permutation is same-as or reverse-of current perm,
      if (any([all(new_permute==permute) all(new_permute==permute(N:-1:1))]))
        diff_obj = init_temp;       % Call it a large difference
        new_obj = obj;
      else
        new_obj = corrobj(r(new_permute,new_permute),w);
        diff_obj = obj - new_obj;   % Else use actual difference
        if (total_perms <= distrib_length)  % Save in distribution
          distrib(total_perms) = new_obj;
        end;
      end;

      % Metropolis algorithm to accept or not
      prob = exp(-diff_obj/temp);
      if (prob > rand)              % If oracle says yes, do it
        permute = new_permute;
        obj = new_obj;
        nsuccess = nsuccess + 1;
        best = addlist(best,[obj permute],0); % Pocket the best so far
      end;

      if (nsuccess >= max_success)  % Finish early if have enough successful
        break;                      %   changes
      end;
    end; % for p=1:max_pertemp

    disp(sprintf( ...
      '  Step=%3.0f, obj=%7.4f, best=%7.4f, replacements=%3.0f, perms=%7.0f', ...
      t, obj, best(1,1), nsuccess, total_perms));
    if (get_stats)
      stats = [stats; t obj best(1,1) total_perms];
    end;

    % If no improvement and have converged on best-so-far, we're done
    if (all([nsuccess==0, permute==best(1,2:N+1)]))
      break;
    end;

    temp = temp * temp_reduce;      % Else reduce temp for next cycle
  end; % for t=1:max_tsteps;

  obj = best(:,1);                  % Vector of best obj-fn values
  permute = best(:,2:N+1);          % Matrix of best permutations
  p = permute(1,:);                 % Best permutation
  rc = r(p,p);                      % Best matrix

  for i=2:nbest                     % Chk congruency of best permutations
    forward = 0;
    for j=1:N                       % Chk pos of each element against
      k = find(p == permute(i,j));  %   %   pos in best permutation
      forward = forward + (j-k)^2;  % Accum differences in position
    end;
    backward = 0;
    for jj=N:-1:1                   % Repeat for reverse permutation
      j = N-jj+1;
      k = find(p == permute(i,jj));
      backward = backward + (j-k)^2;
    end;
    if (backward < forward)         % Reverse permutation if necessary
      permute(i,1:N) = permute(i,N:-1:1);
    end;
  end;

  plot(stats(:,1),stats(:,2),,'k',stats(:,1),stats(:,3),'k');  % Plot convergence
  return;

