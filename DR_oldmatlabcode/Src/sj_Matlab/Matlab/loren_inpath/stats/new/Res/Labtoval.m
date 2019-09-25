% LABTOVAL: Converts a string matrix of labels to a vector of corresponding 
%           numeric values, one value per row.  
%
%           See valtolab().
%
%     Usage: [vals,uvals,ords,ovals] = labtoval(s)
%
%           s = [r x c] string matrix.
%           ----------------------------------------------
%           vals =  [r x 1] vector of corresponding values.
%           uvals = [u x 1] vector of unique sorted values.
%           ords =  [u x 1] matching vector of ordinals [1; 2; ... ].
%           ovals = [r x 1] matching vector of ordinals corresponding to 'vals'.
%

% RE Strauss, 6/10/98

function [vals,uvals,ords,ovals] = labtoval(s)

  n = abs(s);                         % Ascii values of lower-case letters
  [r,c] = size(s);

  z = abs(' ')-1;                     % Effective zero
  b = abs('z')-z+1;                   % Base for numeric conversion

  vals = zeros(r,1);                  % Numeric values for subspp
  for i = 1:c
    vals = vals + (n(:,i)-z) * b.^(c-i);
  end;

  uvals = uniquef(vals,1);
  ords = [1:length(uvals)]';

  ovals = zeros(r,1);
  for i = 1:length(uvals)
    j = find(vals == uvals(i));
    ovals(j) = i * ones(length(j),1);
  end;

  return;

