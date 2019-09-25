% TRIRCRND: Get a closed polygon and generate an indefinite series of random
%           triangulations of it via recursion.
%           Triangulations are written to an output file.
%

% RE Strauss, 10/24/94

function trircrnd
  disp('  Generate random triangulations:');
  N = input('    Number of boundary points: ');
  if (isempty(N))
    return;
  end;
  iterations = input('    Number of iterations: ');
  if (isempty(iterations))
    return;
  end;

  fid = fopen('e:\matlab\output\trircrnd.out','w');
  display_every = 50;

  Ttot = 0;                     % Total triangles generated

  % --- Assign circular coordinates for polygon --- %

  pt = zeros(N,2);
  theta = 2*pi/N;
  for i=1:N
    pt(i,1) = cos((i-1)*theta + pi);
    pt(i,2) = sin((i-1)*theta + pi);
  end;

  % --- Random triangulations --- %

  while (Ttot < iterations)




    fprintf(fid, ' %1.0f', newv);     % Write it to file
    fprintf(fid, '\r\n');

    Ttot = Ttot+1;                    % Increment total produced
    if (any([v ~= newv]))             % Increment total modified
      Tmod = Tmod+1;
    end;

    if (~rem(Ttot,display_every))     % Display every so-many iterations
      Tperc = Tmod * 100/Ttot;
      fprintf('  Total=%1.0f, modified=%1.0f, percent modified=%4.2f\r\n', ...
            Ttot, Tmod, Tperc);
    end;
  end;  % One random triangulation

  fprintf('\r\n  Total=%1.0f\r\n', Ttot);
  fprintf('  Number of points = %1.0f\r\n', N);
  fclose(fid);
  return;




