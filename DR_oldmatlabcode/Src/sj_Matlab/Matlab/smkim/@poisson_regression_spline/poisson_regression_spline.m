classdef poisson_regression_spline
% Two-dimensional tensor product cubic B-spline surface, with periodic boundary
% condition along the 2nd dimension

  properties (GetAccess = public, SetAccess = private)
    grid_size
    knots
    data
    n_data
    n_knots
    % precomputed coefficients
    A
    B
    C
    D
  end

  properties (Constant = true)
  end

  methods (Access = public)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructor
    function self = poisson_regression_spline(knots,data)
      if ~isstruct(knots) || ~isscalar(knots) || ...
          ~all(isfield(knots,{'x','phi'})) || ...
          ~isvector(knots.x) || ~isreal(knots.x) || ...
          ~isfloat(knots.x) || ~all(isfinite(knots.x)) || ...
          (numel(knots.x) < 4) || ~all(diff(knots.x) > 0) || ...
          ~all(abs(diff(knots.x,2)) < sqrt(max(eps(knots.x)))) || ...
          ~isvector(knots.phi) || ~isreal(knots.phi) || ...
          ~isfloat(knots.phi) || ~all(isfinite(knots.phi)) || ...
          (numel(knots.phi) < 4) || ...
          any((knots.phi < -pi) | (knots.phi > +pi)) || ...
          ~all(diff(knots.phi) > 0) || ...
          ~all(abs(diff([diff(knots.phi(:)); 
          knots.phi(1)+2*pi-knots.phi(end)])) < ...
          sqrt(max(eps(knots.phi))))
        error('KNOTS is not valid');
      end
      % Generate two-dimensional grid of knots
      [self.knots.x, self.knots.phi] = ndgrid(knots.x,knots.phi);
      self.grid_size = [numel(knots.x), numel(knots.phi)];
      self.n_knots = prod(self.grid_size);
      % We rely on these properties throughout the methods of this class, so
      % it's good to check that they hold!
      assert(isequal(size(self.knots.x),size(self.knots.phi), ...
          self.grid_size) && ...
          (numel(self.knots.x) == ...
          nnz(bsxfun(@eq,self.knots.x(:,1),self.knots.x))) && ...
          (numel(self.knots.phi) == ...
          nnz(bsxfun(@eq,self.knots.phi(1,:),self.knots.phi))) );

      % Check that data.x, data.phi and data.count are valid column vectors of
      % the same size, that data.timestep is a positive scalar, and that all
      % elements of data.x lie between interior points of knots.x
      if ~isstruct(data) || ~isscalar(data) || ...
          ~all(isfield(data,{'x','phi','count','timestep'})) || ...
          ~isequal(size(data.x,1),size(data.phi,1),size(data.count,1), ...
          numel(data.x),numel(data.phi),numel(data.count)) || ...
          ~isvector(data.x) || ~isreal(data.x) || ...
          ~isfloat(data.x) || ~all(isfinite(data.x)) || ...
          any((data.x < knots.x(2)) | (data.x >= knots.x(end-1))) || ...
          ~isvector(data.phi) || ~isreal(data.phi) || ...
          ~isfloat(data.phi) || ~all(isfinite(data.phi)) || ...
          any((data.phi < -pi) | (data.phi > +pi)) || ...
          ~isvector(data.count) || ~isreal(data.count) || ...
          ~isfloat(data.count) || ~all(isfinite(data.count)) || ...
          ~all(data.count >= 0) || ~all(round(data.count) == data.count) || ...
          ~isscalar(data.timestep) || ~isreal(data.timestep) || ...
          ~isfloat(data.timestep) || ~isfinite(data.timestep) || ...
          ~(data.timestep > 0)
        error('DATA is not valid');
      end
      self.data.x = data.x;
      self.data.phi = data.phi;
      self.data.count = data.count;
      self.data.timestep = data.timestep;
      self.n_data = numel(self.data.count);
      % For each data point, construct lookup table of bounding knots for fast
      % computation of log likelihood, and for each knot, construct lookup table
      % of adjacent knots for fast computation of penalty
      [self.A, self.B, self.C, self.D] = self.precompute();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate of Poisson process intensity at data points, as a function of
    % B-spline coefficients
    rate = lambda(self,coeff);

    % Evaluate negative log-likelihood of B-spline knot coefficients given the
    % observed Poisson counts
    [NLL, NLL_grad, NLL_hess] = neg_log_like(self,coeff);

    % Evaluate rate at knots
    rate = evaluate(self,coeff);

    % Compute sum of squared 2nd-order differences between interior knots along
    % each dimension; these can be used as an approximate curvature penalty. The
    % return value is a two-element row vector; the first element is the penalty
    % term along the x dimension, and the second term is the penalty along the
    % phi dimension
    [P, P_grad] = penalty(self,coeff);

  end

  methods (Access = private)
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each data point, look up nearby knots
    ind = find_adjacent_knots(self);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each data point, look up nearby knots
    [knot_ind, u, v] = construct_lookup_table(self);
  end

  % Set methods
  methods


  end
  

end



