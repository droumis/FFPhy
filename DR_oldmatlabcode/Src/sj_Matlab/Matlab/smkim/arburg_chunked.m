function varargout = arburg_chunked(x, p)
%ARBURG_CHUNKED AR parameter estimation via Burg method, for signal with indeterminate gaps.
%
%   This function is a trivial modification of the ARBURG function included in
%   the MATLAB Signal Processing Toolbox.
%
%   A = ARBURG_CHUNKED(X,ORDER) returns the polynomial A corresponding to the AR
%   parametric signal model estimate of signal X using Burg's method. X is a
%   cell array of vectors, such that X{1}, X{2}, X{3}, ... are non-contiguous
%   segments of the same overall stationary autoregressive process. ORDER is the
%   model order of the AR system.
%
%   [A,E] = ARBURG(...) returns the final prediction error E (the variance
%   estimate of the white noise input to the AR model).
%
%   [A,E,K] = ARBURG(...) returns the vector K of reflection coefficients
%   (parcor coefficients).
%
%   See also ARBURG in the MATLAB Signal Processing toolbox.
%
% Written by SMK, 2009 June 18.
%

error(nargchk(2,2,nargin))

if isempty(p) || ~(p == round(p)) || ~(p > 0)
  error('Model order must be an integer.')
end

if ~iscell(x) || ~all(cellfun(@isvector,x))
  error('X must be a cell array of vectors');
end
for i = 1:numel(x)
  if (length(x{i}) <= 2*p)
    error('Each cell of X must have length greater than 2*ORDER');
  elseif issparse(x{i}),
    error('Input signal cannot be sparse.')
  end
end

% Coerce cells to be column vectors
x = cellfun(@(c) c(:),x,'UniformOutput',false);
% This vector of lengths help us to keep track of the gaps between the segments
% of the signal
n = cellfun(@numel,x);
n = n(:);
% Now concatenate into a single long column vector (we rely on n to recover the
% original cell array)
x = vertcat(x{:});
% Initialization
ef = x;
eb = x;
a = 1;

% Initial error
E = x' * x / sum(n);

% Preallocate 'k' for speed.
k = zeros(1, p);

for m = 1:p
  % These are logical vectors for selecting the desired elements of ef and eb.
  % With care, we achieve logical indexing that respects the cell array
  % structure of the original input
  idx = n - m + 1;
  efp_bitmask = true(size(ef));
  efp_bitmask(1 + [0; cumsum(idx(1:end-1))]) = false;
  ebp_bitmask = true(size(eb));
  ebp_bitmask(cumsum(idx(1:end))) = false;

  % Calculate the next order reflection (parcor) coefficient
  efp = ef(efp_bitmask);
  ebp = eb(ebp_bitmask);
  num = -2 .* ebp' * efp;
  den = efp' * efp + ebp' * ebp;
 
  k(m) = num ./ den;
 
  % Update the forward and backward prediction errors
  ef = efp + k(m) * ebp;
  eb = ebp + k(m)' * efp;
 
  % Update the AR coeff.
  a= [a; 0] + k(m) * [0; conj(flipud(a))];
 
  % Update the prediction error
  E(m+1) = (1 - k(m)' * k(m)) * E(m);
end

a = a(:).'; % By convention all polynomials are row vectors
varargout{1} = a;
if nargout >= 2
    varargout{2} = E(end);
end
if nargout >= 3
    varargout{3} = k(:);
end


% For your reference, this is the original code for the arburg.m from the MATLAB
% Signal Processing Toolbox.
%   Author(s): D. Orofino and R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.12.4.2 $  $Date: 2004/12/26 22:15:20 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%error(nargchk(2,2,nargin))
%
%[mx,nx] = size(x);
%if isempty(x) || min(mx,nx) > 1,
%   error('X must be a vector with length greater than twice the model order.');
%elseif isempty(p) || ~(p == round(p))
%   error('Model order must be an integer.')
%end
%if issparse(x),
%   error('Input signal cannot be sparse.')
%end
%
%x  = x(:);
%N  = length(x);
%
%% Initialization
%ef = x;
%eb = x;
%a = 1;
%
%% Initial error
%E = x'*x./N;
%
%% Preallocate 'k' for speed.
%k = zeros(1, p);
%
%for m=1:p
%   % Calculate the next order reflection (parcor) coefficient
%   efp = ef(2:end);
%   ebp = eb(1:end-1);
%   num = -2.*ebp'*efp;
%   den = efp'*efp+ebp'*ebp;
%   
%   k(m) = num ./ den;
%   
%   % Update the forward and backward prediction errors
%   ef = efp + k(m)*ebp;
%   eb = ebp + k(m)'*efp;
%   
%   % Update the AR coeff.
%   a=[a;0] + k(m)*[0;conj(flipud(a))];
%   
%   % Update the prediction error
%   E(m+1) = (1 - k(m)'*k(m))*E(m);
%end
%
%a = a(:).'; % By convention all polynomials are row vectors
%varargout{1} = a;
%if nargout >= 2
%    varargout{2} = E(end);
%end
%if nargout >= 3
%    varargout{3} = k(:);
%end

