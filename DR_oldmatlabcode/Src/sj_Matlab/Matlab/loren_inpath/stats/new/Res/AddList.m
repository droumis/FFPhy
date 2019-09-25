% ADDLIST:  For optimization problems, adds a row item to an accumulated
%           list of optimal configurations.  Tailored for CORRPERM.
%
%     Usage: [newlist,update] = addlist(list,item,{minflag})
%
%           list =    [n x m] matrix comprising a list of n objects, the first
%                       element of which is the value of the objective function
%                       and the remaining elements the parameters of the
%                       configuration being optimized
%           item =    [1 x m] item to be inserted into the list
%           minflag = boolean flag: (default TRUE)
%                       TRUE (1) =  objective function is being minimized
%                       FALSE (0) = objective function is being maximized
%           newlist = updated [n x m] matrix
%           update =  boolean flag, TRUE if the list was updated, else FALSE
%

% RE Strauss, 3/11/95

function [newlist,update] = addlist(list,item,minflag)
  TRUE = 1; FALSE = 0;

  [N,M] = size(list);                 % Size of list
  [n,m] = size(item);
%  if (n~=1)
%    error('  ADDLIST: item must be row vector');
%  end;
%  if (m~=M)
%    disp([m M]);
%    error('  ADDLIST: item wrong size for list');
%  end;

  newlist = list;
  it = item;
  update = FALSE;

%  if (nargin < 3)                     % Default minflag
%    minflag = 1;
%  end;

%  if (minflag)                        % If minimizing rather than maximizing,
%    list(1,:) = -1*list(1,:);         % make objective-function values negative
%    it(1) = -it(1);
%  end;

  if (it(1,1) <= list(N,1))           % If item not good enough for list,
    return;                           %   forget it
  end;

  % Scan list bottom to top
  for i=N-1:-1:1
    % If obj-fn value already in list, if config or reverse config in list,
    %   forget it
    if (list(i,1)-it(1,1) < 1e-6)
      if (all(list(i,3:M)==it(1,3:M)))
        return;
      elseif (all(list(i,3:M)==it(1,M:-1:3)))
        return;
      end;

    % Find place for new item
    elseif (list(i,1) > it(1,1));
      newlist(i+2:N,:) = newlist(i+1:N-1,:);  % Shift rows behind it
      newlist(i+1,:) = item(1,:);     % Insert item
      update = TRUE;
      return;
    end;
  end;

  newlist(2:N,:) = newlist(1:N-1,:);  % If haven't inserted item,
  newlist(1,:) = item(1,:);           %   insert at front
  update = TRUE;

  return;

