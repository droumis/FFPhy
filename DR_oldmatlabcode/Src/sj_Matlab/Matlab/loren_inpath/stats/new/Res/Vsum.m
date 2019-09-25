% VSUM: sum of a vector, ignoring non-finite values.
%
%     Usage: s = vsum(x)
%
%         x = any vector.
%         -------------------------
%         s = sum of finite values.
%

function s = vsum(x)
  
  s = 0;
  for i = 1:length(x)
    if (isfinite(x(i)))
      s = s + x(i);
    end;
  end;

  return;
  
