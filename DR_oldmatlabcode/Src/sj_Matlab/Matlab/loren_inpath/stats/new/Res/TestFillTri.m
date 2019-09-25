n = 5;
m = zeros(n,n);
counter = 0;

for r = 1:(n-1)
  for c = (r+1):n
    counter = counter + 1;
    m(r,c) = counter;
  end;
end;

% m = m + m';