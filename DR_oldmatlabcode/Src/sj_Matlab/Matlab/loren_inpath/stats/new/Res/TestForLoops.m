
x = zeros(5,5);

counter = 0;

for r = 1:2:5
  for c = 1:5
    counter = counter + 1;
    x(r,c) = counter;
  end;
end;

x


