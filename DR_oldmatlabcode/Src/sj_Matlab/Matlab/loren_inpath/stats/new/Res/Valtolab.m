% VALTOLAB: Converts a vector of numeric values, previously created by function 
%           labtoval(), back to a string matrix of character labels (by row).
%
%           See labtoval().
%
%     Usage: labels = valtolab(vals,lablen)
%
%           vals = vector of numeric values.
%           lablen = length of corresponding character labels.
%           --------------------------------------------------
%           labels = corresponding character-string labels.
%

% RE Strauss, 6/10/98

function labels = valtolab(vals,lablen)

  r = length(vals);
  z = abs(' ')-1;                       % Effective zero
  b = abs('z')-z+1;                     % Base for numeric conversion

  labels = zeros(r,lablen);             % Allocate output string matrix

  uv = uniquef(vals);                   % Unique numeric values
  ulab = zeros(length(uv),lablen);      % Unique numeric values for subspp

  for i = 1:length(uv)                  % Cycle thru unique values
    u = uv(i);
    for j = 1:lablen-1                     
      v = 0;
      while (u > b.^(lablen-j))
        u = u - b.^(lablen-j);
        v = v+1;
      end;
      ulab(i,j) = v + z;
    end;
    ulab(i,lablen) = u + z;

    j = find(vals == uv(i));
    labels(j,:) = ones(length(j),1) * ulab(i,:);
  end;

  labels = char(labels);

  return;

