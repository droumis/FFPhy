function rate = lambda(self,coeff)
%LAMBDA Method for poisson_regression_spline class

  %{
  assert(isequal(numel(coeff),prod(self.grid_size)) && isreal(coeff) && ...
      all(isfinite(coeff(:))));
  %}

  % Given coefficients, compute estimated Poisson rate for each of the data bins
  rate = exp(self.B * coeff);  

  %{
  assert(isequal(numel(rate),size(rate,1),self.n_data));
  %}

end

