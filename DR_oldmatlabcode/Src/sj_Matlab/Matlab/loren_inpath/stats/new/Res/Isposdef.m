% ISPOSDEF: Boolean function that determines whether a square symmetric matrix 
%           is positive-definite (all eigenvalues > zero).  If the matrix is not 
%           square-symmetric, returns 'false' and displays a warning message.
%           See posdef() to find a corresponding positive-definite matrix if 
%           this one is not such.
%
%     Usage: b = isposdef(S)
%
%         S = square symmetric matrix.
%         ----------------------------------
%         b = boolean response (true/false).
%

% RE Strauss, 11/2/01

function b = isposdef(S)
  b = 0;
  if (issqsym(S))
    [evect,eval] = eig(S);                  % Eigen decomposition
    eval = diag(eval);
    if (all(eval>0))                        % Matrix is pos def
      b = 1;
    end;
  else
    disp('  ISPOSDEF warning: matrix is not square-symmetric.');
  end;

  return;
