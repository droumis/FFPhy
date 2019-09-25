function [h, levels] = contourc_lines(c,varargin)
%CONTOURC_LINES Draws line objects corresponding the contour matrix output of CONTOURC.
%
%   [H, LEVELS] = CONTOURC_LINES(C) takes a contour matrix C of the type
%   returned by CONTOURC (see MATLAB documentation) and draws the contours as
%   line objects. H is a vector of handles to the line objects. LEVELS is a
%   vector of the contour levels corresponding to these lines.
%
%

if ~isreal(c) || isempty(c) || (ndims(c) ~= 2) || (size(c,1) ~= 2)
  error('C must be a real 2-by-N matrix');
end

levels = [];
h = [];
i = 1;
while (i < size(c,2))
  levels(end+1) = c(1,i);
  n = c(2,i);
  if ~isfinite(n) || (round(n) ~= n)
    error('C is not a valid contour matrix');
  end
  j = i + (1:n);
  try
    h(end+1) = line(c(1,j),c(2,j),varargin{:});
  catch
    error('Error calling LINE. Invalid arguments?');
  end
  i = i + n + 1;
end


