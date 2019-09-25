% TRAJECT2: Fits two-stage splined ontogenetic trajectories to juveniles
%     and adults, and produces goodness-of-fit statistics.
%
%     Usage: s = traject2(stage,bodysize,char)
%
%        stage =    vector containing group labels for juveniles and adults.
%        bodysize = vector of body sizes for all individuals.
%        char =     [n x p] matrix of characters to be regressed on size.
%        s =        [p x 2] matrix of coordinates of inflection points.
%
%     Written to 'c:\matlab\output\trajcrd1.dat':
%        Coordinates of inflection point, and predicted low & high points
%        for bilinear-spline trajectory.
%     Written to 'c:\matlab\output\trajcrd2.dat':
%        Coordinates of predicted low & high points for pooled linear
%        trajectory.
%     Written to 'c:\matlab\output\trajcoef.dat':
%        Slope and intercept coefficients for bilinear spline and pooled
%        linear regression.
%     Written to 'c:\matlab\output\trajstat.dat':
%        SSEs, improvement-of-fit F-statistic, df's, probability.
%

% RE Strauss, 4/30/94

function s = traject2(stage,bodysize,char)
   iter = input('Number of bootstrap iterations (0=exact solution only): ');
   if (iter>0)                      % Initiate random-number generator
      clk = clock;                  %   for bootstrapping
      rand('seed',clk(6));
   end;

   [j_size,a_size] = splitgrp(stage,bodysize); % Split size-vector by stage
   [j_char,a_char] = splitgrp(stage,char);     % Split char-matrix by stage

   % Allocate vectors and matrices
   s_opt = zeros(1,2);              % Coord vectors
   s0 = zeros(1,2);
   minpred = zeros(1,2);
   maxpred = zeros(1,2);
   minlinear = zeros(1,2);
   maxlinear = zeros(1,2);
   centr_j = zeros(1,2);
   centr_a = zeros(1,2);
   minsize = ones(1,2);             % Augmented vectors
   maxsize = ones(1,2);
   s = [];                          % Accumulation matrices
   b = [];

   % Initializations
   centr_j(1) = mean(j_size);       % Size centroids
   centr_a(1) = mean(a_size);
   min_j = min(j_size);
   max_a = max(a_size);
   extension = 0.05*(max_a - min_j);
   minsize(1) = min_j - extension;  % Extend axes beyond data
   maxsize(1) = max_a + extension;
   minpred(1) = minsize(1);         % Set X-coords of min,max
   maxpred(1) = maxsize(1);
   minlinear(1) = minsize(1);
   maxlinear(1) = maxsize(1);

   % Optimization options
   options = foptions;              % Get default options
   options(9) = 0;                  % Suppress analytic gradient check
   options(14) = 400;               % Max iterations
%  options(1) = 1;                  % Display tabular results

   % Overwrite files: delete if exist
   if (exist('c:\matlab\output\trajcrd1.dat'))
      delete c:\matlab\output\trajcrd1.dat;
   end;
   if (exist('c:\matlab\output\trajcrd2.dat'))
      delete c:\matlab\output\trajcrd2.dat;
   end;
   if (exist('c:\matlab\output\trajcoef.dat'))
      delete c:\matlab\output\trajcoef.dat;
   end;
   if (exist('c:\matlab\output\trajstat.dat'))
      delete c:\matlab\output\trajstat.dat;
   end;

   % Write headings to files
   fprintf('c:\matlab\output\trajcrd1.dat', ...
      '                          Bilinear Spline\n' );
   fprintf('c:\matlab\output\trajcrd1.dat', ...
      '      Inflection Point   Minimum Predicted  Maximum Predicted\n\n' );

   fprintf('c:\matlab\output\trajcrd2.dat', ...
      '          Simple Linear Regression \n' );
   fprintf('c:\matlab\output\trajcrd2.dat', ...
      '    Minimum Predicted  Maximum Predicted\n\n' );

   fprintf('c:\matlab\output\trajcoef.dat', ...
      '        Juv Trajectory   Adult Trajectory   Pooled Trajectory\n' );
   fprintf('c:\matlab\output\trajcoef.dat', ...
      '       Slope Intercept    Slope Intercept    Slope Intercept\n\n' );

   fprintf('c:\matlab\output\trajstat.dat', ...
      '      SSE(F)   SSE(R)      F      df_F df_R   Prob\n\n' );

   [npts,p] = size(char);           % Get numbers of obs and characters
   df_f = npts-4;                   % df for full model (bilinear spline)
   df_r = npts-2;                   % df for reduced model (pooled regr)

   % Cycle through characters
   for c = 1:p
      disp_str = ['   Fitting char ',int2str(c),'...'];
      disp(disp_str);

      J = [j_size j_char(:,c)];     % Compose crd matrices for juvs, adults
      A = [a_size a_char(:,c)];

      % Reduced model: simple linear regression

      [b_lin,stats] = linregr(bodysize,char(:,c)); % Regr pooled data
      r2 = stats(1);
      mse_r = stats(2);
      df_r = stats(4);
      sse_r = mse_r * df_r;

      minlinear(2) = minsize * b_lin;       % Pred min & max values
      maxlinear(2) = maxsize * b_lin;

      % Full model: bilinear spline, seeded from reduced model but using
      % centroid of group centroids rather than grand centroid

      s0 = mean([mean(J); mean(A)]);
      s0(2) = [s0(1) 1] * b_lin';
      p = [s0 b_lin(1) b_lin(1)];                % Seed parameter vector
      p_opt = trajopt('trajfun2',p,options,J,A); % Optimize solution

      s_opt = p_opt(1:2);                    % Retrieve inflection pt
      [sse_f,b_opt] = trajfun2(p_opt,J,A);   % Get stats for fitting
      minpred(2) = minsize * b_opt(1:2)';    % Pred min & max values
      maxpred(2) = maxsize * b_opt(3:4)';

      if (iter > 0)                 % Bootstrapped confidence interval
         disp('     Bootstrapping...');

         Bs = [];                      % Allocate accumulation matrices
         Bb = [];
         Br = [];

         for bit = 1:iter
            BJ = bootsamp(J);
            BA = bootsamp(A);
            JA = [BJ; BA];

            B_lin = linregr(JA(:,1),JA(:,2));  % Simple regr of pooled data

            s0 = mean([mean(BJ); mean(BA)]);   % Bilinear regression
            s0(2) = [s0(1) 1] * B_lin';
            p = [s0 B_lin(1) B_lin(1)];        % Seed parameter vector
            p_opt = trajopt('trajfun2',p,options,BJ,BA); % Optimize solution
            Bs_opt = p_opt(1:2);               % Retrieve inflection pt
            [Bsse_f,Bb_opt] = trajfun2(p_opt,J,A);   % Get stats for fitting

            Bs = [Bs; Bs_opt];         % Accumulate bootstrapped results
            Bb = [Bb; Bb_opt];
            Br = [Br; B_lin];
         end;

         mean_s = mean(Bs);        % Get statistics of bootstrapped
         median_s = median(Bs);    %   distributions
         std_s = std(Bs);

         mean_b = mean(Bb);
         median_b = median(Bb);
         std_b = std(Bb);

         mean_r = mean(Br);
         median_r = median(Br);
         std_r = std(Br);
      end;

      s = [s; s_opt];                  % Accumulate inflections
      b = [b; b_opt];                  %   and regression coefficients

      % Improvement-of-fit F-statistic and corresponding probability
      F = ((sse_r-sse_f)/(df_r-df_f)) / (sse_f/df_f);
      prob = fcdf(F,2,df_f);

      % Write results to files
      fprintf('c:\matlab\output\trajcrd1.dat', ...
         '%3.0f  %8.3f %8.3f  ', c, s_opt(1),s_opt(2) );
      fprintf('c:\matlab\output\trajcrd1.dat', ...
         '%8.3f %8.3f  ', minpred(1),minpred(2) );
      fprintf('c:\matlab\output\trajcrd1.dat', ...
         '%8.3f %8.3f\n', maxpred(1),maxpred(2) );

      if (iter > 0)                 % Bootstrap results
         fprintf('c:\matlab\output\trajcrd1.dat', ...
            'mean %8.3f %8.3f\n', mean_s(1),mean_s(2) );
         fprintf('c:\matlab\output\trajcrd1.dat', ...
            'med  %8.3f %8.3f\n', median_s(1),median_s(2) );
         fprintf('c:\matlab\output\trajcrd1.dat', ...
            'stde %8.3f %8.3f\n\n', std_s(1),std_s(2) );
      end;

      fprintf('c:\matlab\output\trajcrd2.dat', ...
         '%2.0f  %8.3f %8.3f  ', c, minlinear(1),minlinear(2) );
      fprintf('c:\matlab\output\trajcrd2.dat', ...
         '%8.3f %8.3f\n', maxlinear(1),maxlinear(2) );

      fprintf('c:\matlab\output\trajcoef.dat', ...
         '%3.0f  %8.4f %8.4f  ', c, b_opt(1),b_opt(2) );
      fprintf('c:\matlab\output\trajcoef.dat', ...
         '%8.4f %8.4f  ', b_opt(3),b_opt(4) );
      fprintf('c:\matlab\output\trajcoef.dat', ...
         '%8.4f %8.4f\n', b_lin(1),b_lin(2) );

      if (iter > 0)                 % Bootstrap results
         fprintf('c:\matlab\output\trajcoef.dat', ...
            'mean %8.4f %8.4f  ', mean_b(1),mean_b(2) );
         fprintf('c:\matlab\output\trajcoef.dat', ...
            '%8.4f %8.4f  ', mean_b(3),mean_b(4) );
         fprintf('c:\matlab\output\trajcoef.dat', ...
            '%8.4f %8.4f\n', mean_r(1),mean_r(2) );

         fprintf('c:\matlab\output\trajcoef.dat', ...
            'med  %8.4f %8.4f  ', median_b(1),median_b(2) );
         fprintf('c:\matlab\output\trajcoef.dat', ...
            '%8.4f %8.4f  ', median_b(3),median_b(4) );
         fprintf('c:\matlab\output\trajcoef.dat', ...
            '%8.4f %8.4f\n', median_r(1),median_r(2) );

         fprintf('c:\matlab\output\trajcoef.dat', ...
            'stde %8.4f %8.4f  ', std_b(1),std_b(2) );
         fprintf('c:\matlab\output\trajcoef.dat', ...
            '%8.4f %8.4f  ', std_b(3),std_b(4) );
         fprintf('c:\matlab\output\trajcoef.dat', ...
            '%8.4f %8.4f\n\n', std_r(1),std_r(2) );
      end;

      fprintf('c:\matlab\output\trajstat.dat', ...
         '%2.0f  %8.4f %8.4f  ', c, sse_f, sse_r );
      fprintf('c:\matlab\output\trajstat.dat', ...
         '%8.4f   %3.0f  %3.0f  ', F, df_f, df_r );
      fprintf('c:\matlab\output\trajstat.dat', ...
        '%7.5f\n', prob );

   end; % Cycling through characters
   return;
