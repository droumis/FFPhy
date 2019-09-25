function [P, P_grad] = penalty(self,coeff)
%

  rate = reshape(self.D * coeff, self.grid_size);

  % take 2nd-order differences along each dimension, taking care to account for
  % periodicity along phi dimension, and sum the squares of these differences

  P_x = diff(rate,2,1).^2;
  P(1) = sum(P_x(:));

  P_phi = diff([rate, rate(:,[1 2])],2,2).^2;
  P(2) = sum(P_phi(:));

  if (nargout > 1)
    %P_grad = ...
  end

end

