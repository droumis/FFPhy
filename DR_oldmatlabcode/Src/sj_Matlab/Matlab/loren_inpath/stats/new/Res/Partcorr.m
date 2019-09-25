% PARTCORR: Calculates the partial correlation of r12.3, given r12, r13 & r23,
%           and finds the range of possible values of r12, given r13 & r23.
%
%     Usage: [pc,rc] = partcorr(r12,r13,r23)
%
%           r12,r13,r23 = vectors (length n) of correlations.  r12 may be null 
%                           or contain missing (NaN) values.
%           ------------------------------------------------------------------
%           pc =          vector of partial correlations (r12.3) corresponding 
%                           to r12.
%           rc =          [n x 2] matrix of possible ranges of r12, given r13 
%                           & r23.  Col 1 contains min possible r12, col 2 
%                           contains max.
%           

% RE Strauss, 10/6/99

function [pc,rc] = partcorr(r12,r13,r23)
  get_rc = 0;
  if (nargout>1)
    get_rc = 1;
  end;

  if (min(size(r12))>1 | min(size(r13))>1 | min(size(r23))>1)
    error('PARTCORR: input matrices must be vectors');
  end;

  is_row = 0;                             % Flag indicating that input are
  if (size(r12,1)==1)                     %   row vectors
    is_row = 1;
  end;

  r12 = r12(:);                           % Convert to col vectors
  r13 = r13(:);
  r23 = r23(:);

  vectlen = length(r13);
  pc = [];

  r13r23 = r13.*r23;                      % Product of corrs
  D = sqrt((1-r13.^2).*(1-r23.^2));       % Denominator

  if (~isempty(r12))                      % Partial correlations
    if (any(isfinite(r12)))
      pc = (r12-r13r23)./D;
      i = find(~isfinite(pc));
      if (~isempty(i))
        pc(i) = zeros(length(i),1);
      end;
      i = find(pc<-1 | pc>1);
      if (~isempty(i))
        pc(i) = NaN*ones(length(i),1);
      end;
    end;
  end;

  if (get_rc)                             % Constrained ranges
    rc = zeros(vectlen,2);
    rc(:,1) = r13r23-D;
    rc(:,2) = r13r23+D;
  end;

  if (is_row)                             % Transpose to row if necessary
    pc = pc';
  end;

  return;
