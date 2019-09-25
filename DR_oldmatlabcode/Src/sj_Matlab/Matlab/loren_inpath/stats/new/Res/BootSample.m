% BootSample: Produce a bootstrapped sample from a data vector.
%
%     Usage: bx = bootsample(x)
%
%         x = [n x 1] vector.
%         ---------------------------------
%         bx = bootstrapped [n x 1] vector.
%

function bx = bootsample(x)
  n = length(x);
  r = ceil(n*rand(n,1));            % Vector of uniformly distributed integers [1,n]
  bx = x(r);                        % Create bootstrapped sample by sampling x
  
  return;
  