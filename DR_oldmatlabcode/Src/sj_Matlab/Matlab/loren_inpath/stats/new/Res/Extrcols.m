% EXTRCOLS: Extracts up to 20 columns of a matrix into separate column vectors.
%           If the number of output arguments is greater than the number of 
%           input columns, the excess arguments are returned as null matrices.
%
%     Usage: [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,...
%                          v11,v12,v13,v14,v15,v16,v17,v18,v19,v20] = extrcols(x)
%
%           x =       [r x c] matrix.
%           ------------------------------------------------------------------
%           v1-v15 =  [r x 1] column vectors corresponding to first 15 cols of 
%                     input matrix.
%

% RE Strauss, 2/5/00
%   5/30/01 - extended number of output arguments from 12 to 15.

function [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20]...
               = extrcols(x)
  [r,c] = size(x);

  v1 = []; v2 = []; v3 = []; v4 = []; v5 = [];
  v6 = []; v7 = []; v8 = []; v9 = []; v10 = [];
  v11 = []; v12 = []; v13 = []; v14 = []; v15 = [];
  v16 = []; v17 = []; v18 = []; v19 = []; v20 = [];

  if (nargout>=1  & c>=1),  v1  = x(:,1);  end;
  if (nargout>=2  & c>=2),  v2  = x(:,2);  end;
  if (nargout>=3  & c>=3),  v3  = x(:,3);  end;
  if (nargout>=4  & c>=4),  v4  = x(:,4);  end;
  if (nargout>=5  & c>=5),  v5  = x(:,5);  end;
  if (nargout>=6  & c>=6),  v6  = x(:,6);  end;
  if (nargout>=7  & c>=7),  v7  = x(:,7);  end;
  if (nargout>=8  & c>=8),  v8  = x(:,8);  end;
  if (nargout>=9  & c>=9),  v9  = x(:,9);  end;
  if (nargout>=10 & c>=10), v10 = x(:,10); end;
  if (nargout>=11 & c>=11), v11 = x(:,11); end;
  if (nargout>=12 & c>=12), v12 = x(:,12); end;
  if (nargout>=13 & c>=13), v13 = x(:,13); end;
  if (nargout>=14 & c>=14), v14 = x(:,14); end;
  if (nargout>=15 & c>=15), v15 = x(:,15); end;
  if (nargout>=16 & c>=16), v16 = x(:,16); end;
  if (nargout>=17 & c>=17), v17 = x(:,17); end;
  if (nargout>=18 & c>=18), v18 = x(:,18); end;
  if (nargout>=19 & c>=19), v19 = x(:,19); end;
  if (nargout>=20 & c>=20), v20 = x(:,20); end;

  return;
