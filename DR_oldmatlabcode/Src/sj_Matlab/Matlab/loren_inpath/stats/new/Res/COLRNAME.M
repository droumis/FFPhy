% COLRNAME: Given a 3-element RGB color specification, returns either the 
%           standard MATLAB long name of the color or the descriptor 'mixed'.
%
%     Usage: cname = colrname(rgb)
%

% RE Strauss, 5/1/95

function cname = colrname(rgb)
  if (rgb(1)~=1 & rgb(1)~=0)
    cname = 'mixed  ';
    return;
  end;
  if (rgb(2)~=1 & rgb(2)~=0)
    cname = 'mixed  ';
    return;
  end;
  if (rgb(3)~=1 & rgb(3)~=0)
    cname = 'mixed  ';
    return;
  end;

  if (rgb(1)==1)
    if (rgb(2)==1)
      if (rgb(3)==1)
        cname = 'white  ';
      else
        cname = 'yellow ';
      end;
    else
      if (rgb(3)==1)
        cname = 'magneta';
      else
        cname = 'red    ';
      end;
    end;
  else
    if (rgb(2)==1)
      if (rgb(3)==1)
        cname = 'cyan   ';
      else
        cname = 'green  ';
      end;
    else
      if (rgb(3)==1)
        cname = 'blue   ';
      else
        cname = 'black  ';
      end;
    end;
  end;
  return;
