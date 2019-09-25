% FINDMATCH: Given two vectors, reduces them to their common values.
%            If either vector has duplicate values, only the first is kept.
%
%     Usage: [reduced,ix,iy] = findmatch(x,y)
%
%         x,y = vectors of arbitrary length.
%         ----------------------------------------------------------
%         reduced = reduced indentical x and y vectors.
%         ix,iy = corresponding vectors of indices of xr into x, and 
%                 yr into y.
%

% RE Strauss, 9/12/01

function [reduced,ix,iy] = findmatch(x,y)
  [isvectx,ncellsx,iscolx] = isvector(x);
  [isvecty,ncellsy,iscoly] = isvector(y);

  if (~isvectx | ~isvecty)
    error('  FINDMATCH: one or both input matrices not in vector form.');
  end;

  x = x(:);
  y = y(:);

  ix = [1:length(x)]';
  iy = [1:length(y)]';

  [u,f] = uniquef(x);
  if (any(f>1))
    ifreq = find(f>1);
    for i = 1:length(ifreq)
      j = find(x==u(ifreq(i)));
      x(j(2:length(j))) = [];
      ix(j(2:length(j))) = [];
    end;
  end;

  [u,f] = uniquef(y);
  if (any(f>1))
    ifreq = find(f>1);
    for i = 1:length(ifreq)
      j = find(y==u(ifreq(i)));
      y(j(2:length(j))) = [];
      iy(j(2:length(j))) = [];
    end;
  end;

  xlen = length(x);
  ylen = length(y);

  if (xlen < ylen)
    a = x;
    b = y;
    ia = ix;
    ib = iy;
  else
    a = y;
    b = x;
    ia = iy;
    ib = ix;
  end;

  i = isin(a,b);
  iia = find(i>0);
  a = a(iia);
  ia = ia(iia);

  i = isin(b,a);
  iib = find(i>0);
  b = b(iib);
  ib = ib(iib);

  if (xlen < ylen)
    reduced = a;
    ix = ia;
    iy = ib;
  else
    reduced = b;
    ix = ib;
    iy = ia;
  end;

  if (~iscolx)
    reduced = reduced';
    ix = ix';
    iy = iy';
  end;

  return;
