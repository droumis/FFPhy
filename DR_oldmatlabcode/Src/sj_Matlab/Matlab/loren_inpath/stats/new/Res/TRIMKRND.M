% TRIMKRND: Get a closed polygon and generate a series of random
%           triangulations of it via random vertex-count vectors.  Plot them 
%           to the screen.
%

% RE Strauss, 8/24/94

function trimkrnd

  disp('  Digitize a single polygon...');
  [polygon] = getpolyg;
  coords = polygon;
  [N,p] = size(coords);
  coords(N,:) = [];             % Delete closure point
  N = N-1;

  % --- Setup plot for random triangulations --- %

  figure(1);
  clg;

  plot_rows = 4;
  plot_cols = 4;

  % --- Random triangulations via vertex counts --- %

  for nplot = 1:(plot_rows*plot_cols)
%    disp(nplot);

    Bn = 7;                           % Binomial parameters
    Bmean = 2-(6/N);
    Bp = Bmean / Bn;
    v = binornd(Bn,Bp, 1,N)+2;        % Vector of random vertex counts

    sumv = sum(v);                    % Check sum against constraint
    targetsum = 4*N-6;                %   and adjust if necessary

    if (sumv > targetsum)
    diff = sumv-targetsum;            % Difference between sum and target
    index = fix(unifrnd(1,N+1-eps, diff,1));  % Random vertex indices
    while (diff > 0)                  % While difference exists...
        i = index(diff);
        while (v(i) <= 2)             % Find random vertex with count > 2
          i = wrap(i+1,N);
        end;
        v(i) = v(i)-1;                % Decrement vertex count
        diff = diff-1;
      end;
    elseif (sumv < targetsum)
      diff = targetsum-sumv;          % Difference between sum and target
      index = fix(unifrnd(1,N+1-eps, diff,1));  % Random vertex indices
      while (diff > 0)                % While difference exists...
        i = index(diff);
        v(i) = v(i)+1;                % Increment random vertex count
        diff = diff-1;
      end;
    end;

    [v,edge] = trirecon(v,coords);    % Reconstruct triangulation

    % --- Plot random triangulation --- %

    subplot(plot_rows,plot_cols,nplot);
    hold on;

    axis([0 1 0 1]);
    axis('square');
    axis('off');

    for i=1:N
      plot(coords(i,1),coords(i,2),'+r');
    end;

    for i=1:N-1;
      for j=i+1:N
        if (edge(i,j))
          x = [coords(i,1), coords(j,1)];
          y = [coords(i,2), coords(j,2)];
          plot(x,y,'-w');
        end;
      end;
    end;

    merit = -trifunc(edge,[],[1]);    % Add total edge length
    s = sprintf('F=%4.2f', merit);
    title(s);

    hold off;

  end;  % One random triangulation

  return;




