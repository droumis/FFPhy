% CHARSHAPE: Given a measure of general size and scores for a character,
%            finds the projection scores on the minor axis (=PC2), which
%            are a measure of size-independent character shape.  Both the
%            size scores and data matrix are assumed to be log-transformed.
%
%     Usage: Xs = charshape(X,size)
%
%         X =    [n x p] data matrix.
%         size = vector of size scores for each individual.
%         ----------------------------------------------------------------
%         Xs =   corresponding [n x p] matrix of shape scores.
%

% RE Strauss, 9/18/01

function Xs = charshape(X,size)
  [n,p] = size(X);
  [isvect,ncells,iscol] = isvector(X);
  if (isvect & ~iscol)
    X = X';
  end;

  size = size(:);
  ns = length(size);
  if (ns ~= n)
    error('  CHARSHAPE: input matrices not compatible in size.');
  end;
  if (any(~isfinite(size)))
    error('  CHARSHAPE: size vector contains missing values.');
  end;

  Xs = zeros(n,p);                    % Allocate output matrix

  for ip = 1:p                        % For each character,
    x = X(:,ip);                        % Isolate column from data matrix
    s = size;
    j = [1:n]';
    i = find(~isfinite(x));             % Locate any missing values for char
    if (~isempty(i))                    % If have missing values, omit
      x(i) = [];
      s(i) = [];
      j(i) = [];
    end;  

    [loadings,percvar,scores] = pcacov([s,x],2);

    Xs(:,ip) = NaN*ones(n,1);           % Fill in missing values in output matrix
    Xs(j,ip) = scores(:,2);             % Stash non-missing values 
  end;

  return;
