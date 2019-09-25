x = [1 4 6 5 7 9 2 3 2 3 6 7]

for i = length(x):-1:1
  if (x(i) < 5)
    x(i) = 0;
  end;
end;

x
