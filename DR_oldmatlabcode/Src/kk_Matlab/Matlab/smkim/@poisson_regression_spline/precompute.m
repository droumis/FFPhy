function [A, B, C, D] = precompute(self)
%Precompute Method for poisson_regression_spline class

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Part 1: precompute coefficients for computing log-likelihood of knot
  % coefficients given data

  % This is the blending matrix for a cubic B-spline with uniform knots.
  M = [ -1, +3, -3, +1; ...
        +3, -6, +3,  0; ...
        -3,  0, +3,  0; ...
        +1, +4, +1,  0 ] ./ 6;
  % i and j are subscripts into self.knots.x and self.knots.phi
  % Each row of i/j looks up subscripts of the 16 bounding knots for a single
  % data point
  assert(isequal(size(M),[4 4]));
  i = zeros([self.n_data, 16]);
  j = zeros([self.n_data, 16]);
  % For each data point, find the grid of 16 bounding knots
  for n = 1:self.n_data
    % Populate i(n,:), which subscripts along the x dimension
    ii = find(self.data.x(n) >= self.knots.x(:,1),1,'last');
    assert(isscalar(ii) && (ii < self.grid_size(1)-1) && (ii > 1));
    i(n,:) = ii + [ -1, -1, -1, -1, ...
        0, 0, 0, 0, ...
        +1, +1, +1, +1, ...
        +2, +2, +2, +2 ];
    % Populate j(n,:), which subscripts along the phi dimension
    jj = find(self.data.phi(n) >= self.knots.phi(1,:),1,'last');
    if ~isscalar(jj)
      jj = self.grid_size(2);
    end
    j(n,:) = jj + [ -1, 0, +1, +2, ...
        -1, 0, +1, +2, ...
        -1, 0, +1, +2, ...
        -1, 0, +1, +2 ];
  end
  % Enforce periodic boundary conditions along phi dimension
  j = 1+mod(j-1,self.grid_size(2));
  assert(all((j(:) >= 1) & (j(:) <= self.grid_size(2))));
  % For each data point, compute fractional distance between the two nearest
  % bounding knots along each dimension, expand this fractional distance as a
  % cubic polynomial, and multiply by the appropriate blending weights
  u = double( (self.data.x - self.knots.x(i(:,5),1)) ./ ...
      (self.knots.x(i(:,9),1) - self.knots.x(i(:,5),1)) );
  % Tile u so that u(n,:) matches i(n,:)
  u = [u.^3, u.^2, u, ones(size(u))] * M;
  u = [u(:,[1 1 1 1]), u(:,[2 2 2 2]), u(:,[3 3 3 3]), u(:,[4 4 4 4])];
  % We know that the rows of self.knots.phi are monotonically increasing and
  % uniformly spaced, becauase this condition was checked in the constructor
  % call. But let's check the first row again:
  assert( all(diff(self.knots.phi(1,:)) > 0) && ...
      all(abs(diff([diff(self.knots.phi(1,:)), ...
      self.knots.phi(1,1) + 2*pi - self.knots.phi(1,end)])) < ...
      sqrt(max(eps(self.knots.phi(1,:))))) );
  phi_step = mean(diff(self.knots.phi(1,:)));
  % Fractional distance between knots along the phi dimension must wrap around
  % the -pi/+pi boundary
  v = double( mod(self.data.phi - self.knots.phi(1,j(:,2))',phi_step) ./ ...
      mod(self.knots.phi(1,j(:,3)) - self.knots.phi(1,j(:,2)),2*pi)' );
  % Tile v so that v(n,:) matches j(n,:)
  v = repmat([v.^3, v.^2, v, ones(size(v))] * M, [1 4]);
  % Unroll everything to be column vectors (for compatibility with accumarray)
  i = i(:);
  j = j(:);
  u = u(:);
  v = v(:);
  % A is a vector of size [1, self.n_knots]. Given a [self.n_knots, 1] vector of
  % knot coefficients, (A * coeff) is a scalar equal to the first term of the
  % negative log likelihood
  A = accumarray( ...
      sub2ind(self.grid_size,i,j), ...
      -repmat(self.data.count,[16 1]) .* u .* v, ...
      [self.n_knots 1],@sum)';
  % B is a *sparse* matrix of size [self.n_data, self.n_knots], with exactly 16
  % non-zero elements in each row. Given a [self.n_knots, 1] vector of knot
  % coefficients, exp(B * coeff) is a [self.n_data, 1] vector of instantaneous
  % Poisson intensities in each of the self.n_data time bins
  B = accumarray( ...
      { repmat((1:self.n_data)',[16 1]), sub2ind(self.grid_size,i,j) }, ...
      u .* v, ...
      [self.n_data, self.n_knots],@sum,0,true);
  % C is a *sparse* vector of size [1, (self.n_knots).^2], in which each element
  % corresponds to an ordered pair of knots
  C = B' * B;
  C = C(:)';
  clear('i','j','u','v');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Part 2: precompute coefficients for computing rate at each knot location;
  % this is needed for the approximate curvature penalty
    
  % This is the blending matrix for computing the value at a given knot given
  % coefficients in a local 3-by-3 grid
  W = [ 1,  4,  1; ...
        4, 16,  4; ...
        1,  4,  1 ] / 36; 
  [i, j] = ndgrid(1:self.grid_size(1),1:self.grid_size(2));
  i = i(:);
  j = j(:);
  i = [i-1; i-1; i-1; i; i; i; i+1; i+1; i+1];
  j = [j-1; j; j+1; j-1; j; j+1; j-1; j; j+1];
  % Endpoints knots along the x dimension are implicitly duplicated in this
  % lookup
  i(i < 1) = 1;
  i(i > self.grid_size(1)) = self.grid_size(1);
  assert(all( (i(:) >= 1) & (i(:) <= self.grid_size(1)) ));
  % Enforce periodic boundary conditions along phi dimension
  j = 1+mod(j-1,self.grid_size(2));
  assert(all( (j(:) >= 1) & (j(:) <= self.grid_size(2)) ));
  assert(isequal(size(i,1),size(j,1),numel(i),numel(j),9*self.n_knots));

  % D is a *sparse* square matrix of size [self.n_knots, self.n_knots], with
  % exactly 9 non-zero elements in each row. Given a [self.n_knots, 1] vector of
  % knot coefficients, exp(C * coeff) is a [self.n_knots, 1] vector of Poisson
  % intensities at each of the knot locations.
  D = accumarray( ...
      {repmat((1:self.n_knots)',[9 1]), sub2ind(self.grid_size,i,j)}, ...
      kron(W(:),ones([self.n_knots, 1])), ...
      [self.n_knots, self.n_knots],@sum,0,true);

end

