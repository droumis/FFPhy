% MISSNORM: Test the effect of missing-data imputation on the Mahalanobis
%           distance between two groups, as a function of number of characters,
%           degree of correlation among variables, distance between groups,
%           proportion of missing data, and orientation between groups.
%           Data generated are multivariate normal.
%
%           Requires no input.  Saves results in a matrix, and saves session
%           to a MAT file on disk for later recovery.
%

% RE Strauss, 3/10/98
%   11/29/99 - changed calling sequence of mahal().

function results = missnorm
  n = 50;                             % Within-group sample size
  n2 = n*2;
  grp = [ones(n,1); 2*ones(n,1)];
  boot_iter = 100;                    % Bootstrap iterations

  nvar_level = [2 5 10 15 20];        % Numbers of variables
  r_level = [0.1 0.5 0.7 0.9];        % Target correlation levels (same for all vars)
  d2_level = [0:1:5];                 % Assigned Mahalanobis D2 diffs between grps
  theta_level = [0 45 90];            % Orientations: 0 = along PC1, 90 = along PC2
  md_level = [0 0.01 0.02 0.03 0.04 0.06 0.08 0.10];  % Proportions of missing data

  nvar_n = length(nvar_level);
  r_n = length(r_level);
  d2_n = length(d2_level);
  theta_n = length(theta_level);
  md_n = length(md_level);

  prod_n = nvar_n * r_n * d2_n * theta_n * md_n;
  disp(sprintf('Total number of combinations = %6.0f\n', prod_n));

  results = zeros(prod_n,25);   % col 1 = number of variables
                                %     2 = correlation
                                %     3 = variance among randomized corrs
                                %     4 = assigned D2
                                %     5 = proportion of missing data
                                %     6 = theta between group means
                                %     7 = reduced sample size n1
                                %     8 = reduced sample size n2
                                %     9 = D2 for reduced matrix (nonmissing vals only)
                                %    10 =   lower 95% CL (bootstrapped)
                                %    11 =   upper 95% CL
                                %    12 = D2 using PCA algorithm (bootstrapped median)
                                %    13 =   lower 95% CL (bootstrapped)
                                %    14 =   upper 95% CL
                                %    15 =   mean dev for predicted vals
                                %    16 =     stdev dev for predicted vals
                                %    17 =   mean abs dev for predicted vals
                                %    18 =     stdev abs dev for predicted vals
                                %    19 = D2 using EM algorithm (bootstrapped median)
                                %    20 =   lower 95% CL (bootstrapped)
                                %    21 =   upper 95% CL
                                %    22 =   mean dev for predicted vals
                                %    23 =     stdev dev for predicted vals
                                %    24 =   mean abs dev for predicted vals
                                %    25 =     stdev abs dev for predicted vals

  repl = 0;                         % Current combination of parameter values

  for nvar = nvar_level;            % Next value for number of variables
    mu1 = zeros(1,nvar);              % Center 1st group on origin
    mu2 = mu1;                        % Initially center 2nd group on origin
    S = ones(1,nvar);                 % Vector of unit standard deviations

    for r = r_level;                  % Next correlation level
      for assigned_d2 = d2_level;       % Next value of assigned D2 distance
        for md = md_level;                % Next level for proportion missing vals

          R = r*ones(nvar,nvar);            % Correlation matrix
          for i = 1:nvar;                   % Put 1's on diagonal
            R(i,i) = 1;
          end;

          x1 = randmvn(n,mu1,S,R);          % Get data for 1st group: MVN
          x2 = randmvn(n,mu2,S,R);          %   and for 2nd group

          if (nvar>2)
            r1 = corrcoef(x1);              % Correlation matrices of randomized data
            r2 = corrcoef(x2);
            r1 = r1(abs(tril(ones(nvar,nvar))-1));
            r2 = r2(abs(tril(ones(nvar,nvar))-1));
            z1 = 0.5 * log((1+r1)./(1-r1)); % Fisher z-transform of corrs
            z2 = 0.5 * log((1+r2)./(1-r2));
            varz = (var(z1)+var(z2))/2;
          else
            varz = 0;
          end;

          x1 = x1 - ones(n,1)*mean(x1);     % Set 1st group to origin

          X = [x1; x2];                     % Concatenate the two groups

          [evects,evals] = eigen(cov(X));   % PCA of pooled groups
          scores = score(X,evects,nvar);    % Scores for subset of PCs

          s1 = std(scores(:,1));            % Stdev for PC1 scores
          s2 = std(scores(:,2));            %           PC2 scores

          for theta = theta_level;          % Next value of angle between PC1/PC2
            repl = repl+1;                    % Next available row in output matrix
            disp(sprintf('Repl = %5.0f of %5.0f', repl, prod_n));

            x2 = x2 - ones(n,1)*mean(x2);     % Set 2nd group to origin

            th = theta * 2*pi / 360;          % Convert from deg to radians
            a = sqrt(assigned_d2*(s1*s1));    % D2 distance along PC1
            b = sqrt(assigned_d2*(s2*s2));    % D2 distance along PC2
            pc1 = a * cos(th);                % Adjust projections onto PC1,PC2
            pc2 = b* sin(th);                 %   depending on angle theta

            mu2_score = zeros(1,nvar);        % Means of grp 2 on PCs
            mu2_score(1) = pc1;
            mu2_score(2) = pc2;

            mu2 = mu2_score * inv(evects);    % Convert PC means to char means
            x2 = x2 + ones(n,1)*mu2;          % Translate grp 2 to new centroid

            X = [x1; x2];                     % Concatenate the two groups

            Xm = X;                           % Make copy of full data matrix
            nmiss = ceil(md*n2*nvar);         % Total number of missing values
            for i = 1:nmiss                   % Randomly assign missing values
              c = ceil(rand*nvar);            % Choose random column (var)
              plug_value = 1;
              while (plug_value)
                rr = ceil(rand*n2);           % Choose random row (obs)
                if (sum(finite(Xm(rr,:)))>1 & finite(Xm(rr,c))) % If have 2+ good vals
                  Xm(rr,c) = NaN;             %   replace one with missing value
                  plug_value = 0;
                end;
              end;
            end;

            i = finite(Xm);                   % Logical matrix of nonmissing vals
            j = prod(i')';
            Xred = Xm(j,:);                   % Reduced data matrix (good vals only)
            grpred = grp(j);                  % Corresponding group id's

            results(repl,1) = nvar;           % Stash current parameter values
            results(repl,2) = r;
            results(repl,3) = varz;
            results(repl,4) = assigned_d2;
            results(repl,5) = md;
            results(repl,6) = theta;

            N1 = sum(grpred == 1);                                 % Reduced sample sizes
            N2 = sum(grpred == 2);
            results(repl,7) = N1;
            results(repl,8) = N2;

            % Mahalanobis distance, adjusted for number of characters
            if (N1>nvar+1 & N2>nvar+1)                             % If sample sizes sufficient,

              [mahald2,mahald2_ci] = mahal(Xred,grpred,boot_iter);   % Reduced data matrix

              results(repl,9) = mahald2_ci(4,1)/nvar;                       % including only
              results(repl,10:11) = [mahald2_ci(2,1) mahald2_ci(1,2)]/nvar; % complete data
            else
              results(repl,9:11) = [-1 -1 -1];                       % Flag small samp sizes
            end;

            if (N1>nvar+1 & N2>nvar+1)
              if (nmiss>0)
                Xest1 = misspca(Xm(1:n,:));                          % PCA algorithm
                Xest2 = misspca(Xm((n+1):n2,:));                     %   by group
                Xest = [Xest1; Xest2];

                dif = X(:)-Xest(:);
                idif = find(dif~=0);
                dif = dif(idif);

                mean_dev = mean(dif);
                std_dev = std(dif);
                mean_abs_dev = mean(abs(dif));
                std_abs_dev = std(abs(dif));
              else
                Xest = Xm;
                mean_dev = 0;
                std_dev = 0;
                mean_abs_dev = 0;
                std_abs_dev = 0;
              end;

              % Mahalanobis distance, adjusted for number of characters
              [mahald2,mahald2_ci] = mahal(Xest,grp,boot_iter);    % Imputed matrix
              results(repl,12) = mahald2_ci(4,1)/nvar;
              results(repl,13:14) = [mahald2_ci(2,1) mahald2_ci(1,2)]/nvar;

              results(repl,15) = mean_dev;
              results(repl,16) = std_dev;
              results(repl,17) = mean_abs_dev;
              results(repl,18) = std_abs_dev;
            else
              results(repl,12:18) = [-1 -1 -1 -1 -1 -1 -1];        % Flag small samp sizes
            end;

            if (N1>1 & N2>1)
              if (nmiss>0)
                Xest1 = missem(Xm(1:n,:));                             % EM algorithm
                Xest2 = missem(Xm((n+1):n2,:));                        %   by group
                Xest = [Xest1; Xest2];

                dif = X(:)-Xest(:);
                idif = find(dif~=0);
                dif = dif(idif);

                mean_dev = mean(dif);
                std_dev = std(dif);
                mean_abs_dev = mean(abs(dif));
                std_abs_dev = std(abs(dif));
              else
                Xest = Xm;
                mean_dev = 0;
                std_dev = 0;
                mean_abs_dev = 0;
                std_abs_dev = 0;
              end;

              % Mahalanobis distance, adjusted for number of characters
              [mahald2,mahald2_ci] = mahal(Xest,grp,boot_iter);      % Imputed matrix
              results(repl,19) = mahald2_ci(4,1)/nvar;
              results(repl,20:21) = [mahald2_ci(2,1) mahald2_ci(1,2)]/nvar;

              results(repl,22) = mean_dev;
              results(repl,23) = std_dev;
              results(repl,24) = mean_abs_dev;
              results(repl,25) = std_abs_dev;
            else
              results(repl,19:25) = [-1 -1 -1 -1 -1 -1 -1];        % Flag small samp sizes
            end;

            disp(sprintf('  Variables = %2.0f', nvar));            % Write status to screen
            disp(sprintf('  Corr = %4.2f, var(corr) = %6.4f', r, varz));
            disp(sprintf('  Assigned D2 = %3.1f', assigned_d2));
            disp(sprintf('  Missing-data level = %4.2f', md));
            disp(sprintf('  Theta = %2.0f degrees', theta));
            disp(sprintf('  Bootstrapped estimates:'));
            disp(sprintf('    Reduced data:  D2 = %8.4f (%8.4f, %8.4f), N = (%2.0f, %2.0f)', ...
                            results(repl,[9 10 11 7 8])));
            disp(sprintf('    PCA algorithm: D2 = %8.4f (%8.4f, %8.4f), dev = %8.4f', ...
                            results(repl,[12 13 14 17])));
            disp(sprintf('    EM algorithm:  D2 = %8.4f (%8.4f, %8.4f), dev = %8.4f\n', ...
                            results(repl,[19 20 21 24])));

          end; % for theta
          save missnorm;                        % Save current results to disk

        end; % for md_level
      end; % for assigned_d2
    end; % for r
  end; % for nvar

  return;
