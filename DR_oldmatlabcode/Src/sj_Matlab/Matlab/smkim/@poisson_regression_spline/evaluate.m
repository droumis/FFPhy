function rate = evaluate(self,coeff)
% Evaluate Poisson intensity at knots
 
  %{ 
  assert(isequal(numel(coeff),prod(self.grid_size)) && isreal(coeff) && ...
      all(isfinite(coeff(:))));
  %}

  rate = exp(self.D * coeff);

end

