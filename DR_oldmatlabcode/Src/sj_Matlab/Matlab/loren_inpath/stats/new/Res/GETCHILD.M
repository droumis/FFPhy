% GETCHILD: Lists the handle and color of each child of the current plot, 
%           and optionally returns the vector of handles.
%
%     Usage: h = getchild
%

% RE Strauss, 5/1/95

function hout = getchild
  fprintf('\n        Figure     Object   Type     Color\n');
  blanks = '  ';
  hout = [];
  k = 0;

  f = get(gcf,'children');
  for j=1:length(f)
    h = get(f(j),'children');
    hout = [hout; h];
    
    for i=1:length(h)
      t = get(h(i),'type');
      if (t(1:4)=='line')
        c = get(h(i),'color');
        cname = colrname(c);
      elseif (t(1:5)=='patch')
        c = get(h(i),'facecolor');
        cname = colrname(c);
      else
        cname = ' ';
      end;

      k = k+1;
      fprintf(' %3.0f  %f  %f  ', k,f(j),h(i));
      tname = '       ';
      tname(1:length(t)) = t;
      disp([tname blanks cname]);
    end;
  end;
  return;
