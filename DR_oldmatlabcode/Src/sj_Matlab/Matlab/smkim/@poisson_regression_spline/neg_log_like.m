function [NLL, NLL_grad, NLL_hess] = neg_log_like(self,coeff)
%NEG_LOG_LIKE Method for poisson_regression_spline class

  % Given knot coefficients, compute negative log likelihood over all data
  NLL = (self.A * coeff) + self.data.timestep .* sum(exp(self.B * coeff),1); 

  if (nargout > 1)
    % Gradient of negative log likelihood. size(NLL_grad) is [self.n_knots, 1];
    % NLL_grad(i) is the partial derivative of NLL with respect to coeff(i)
    NLL_grad = ( self.A + self.data.timestep .* ...
        sum(bsxfun(@times,self.B,exp(self.B * coeff)),1) )';
  end

  if (nargout > 2)
    % Hessian of negative log likelihood. size(NLL_hess) is 
    % [self.n_knots, self.n_knots]. NLL_hess(i,j) is the second partial
    % derivative of NLL with respect to coeff(i) and coeff(j).
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Computing the Hessian uses much memory!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % This is what we want to do, but unfortunately MATLAB does not yet support
    % sparse arrays of dimension greater than 2
    %{
    NLL_hess = self.data.timestep .* sum(bsxfun(@times, ...
        reshape(exp(self.B * coeff),[1 1 self.n_data]), ...
        self.B' * self.B),3);
    %}
    % This workaround unrolls the sparse matrix as a vector
    NLL_hess = self.data.timestep .* reshape( ...
        sum(bsxfun(@times,exp(self.B * coeff),self.C),1), ...
        [self.n_knots, self.n_knots]);

  end

end

