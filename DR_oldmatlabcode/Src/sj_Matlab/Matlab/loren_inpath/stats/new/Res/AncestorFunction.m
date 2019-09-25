% AncestorFunction: Given either an existing ancestor function of improper form (e.g.,
%         lacking a zero for the root, having the root specified in other than the 
%         last position, or containing polytomies), or a matrix containing 
%         ancestor-descendend-branch length information, returns a proper ancestor 
%         function and other possible vectors containing terminal-taxon labels and 
%         branch lengths.  Fully resolves polytomies using zero branch lengths.
%
%     Usage: [anc_out,brlen,ttlabels] = ancestorfunction(anc,{brlen})
%
%         anc =      an ancestor-function vector, or a 2- or 3-column matrix specifying
%                      descendant labels (col 1), ancestor labels (col 2), and optional
%                      branch lengths (col 3).
%         brlen =    optional corresponding branch-length vector.
%         -----------------------------------------------------------------------------
%         anc =      revised ancestor function in proper form.
%         brlen =    corresponding branch lengths.
%         ttlables = vector of labels for terminal taxa.
%

% RE Strauss, 10/18/02
%   10/21/02 - corrected problem with branch-length vector.

function [anc,brlen,ttlabels] = ancestorfunction(anc,brlen)
  if (nargin<2) brlen = []; end;
  
  patch_anc = 0;
  convert_anc = 0;
  given_brlen = 0;
  ttlabels = [];

  if (~isempty(brlen))
    given_brlen = 1;
  end;
  if (isvector(anc))
    patch_anc = 1;
    if (~given_brlen)
      brlen = ones(size(anc));
    end;
  elseif (size(anc,2)==2)
    convert_anc = 1;
    topo = anc;
    desc = anc(:,2);
    anc = anc(:,1);
    brlen = ones(length(anc));
  elseif (size(anc,2)==3)
    convert_anc = 1;
    desc = anc(:,2);
    brlen = anc(:,3);
    anc = anc(:,1);
    given_brlen = 1;
  else
    error('  AncestorFunction: invalid input matrix.');
  end;
    
  if (convert_anc)
    tt = find(~isin(anc,desc));                     % Find terminal taxa
    nd = find(isin(anc,desc));                      % Find nodes in first column
    ttlabels = anc(tt)';                            % Vector of terminal-taxon labels
    ntt = length(tt);                               % Number of terminal taxa
    anc(tt) = [1:ntt]';                             % Replace tt labels with initial integers
    und = uniquef([anc(nd);desc]);                  % List of unique node labels
    rund = ranks(und) + ntt;                        % Convert to sequential integers following tt labels
    anc(nd) = replace(anc(nd),und,rund);  
    desc = replace(desc,und,rund);  
    anc_new = zeros(1,length(und));                 % Create ancestor function
    brlen_new = zeros(1,length(und));
    for i = 1:length(anc)
      anc_new(anc(i)) = desc(i);
      brlen_new(anc(i)) = brlen(i);
    end;
    if (~any(anc_new)==0)                           % Append a root if none is present
      anc_new = [anc_new, 0];
      brlen_new = [brlen_new, 0];
    end;
    anc = anc_new;
    brlen = brlen_new;
    if (length(und)~=(2*ntt-1))
      patch_anc = 1;
    end;
  end;

  if (patch_anc)                                  % Patch the ancestor function
    if (~any(anc==0))                               % If no zero for root, create one
      anc = [anc, 0];
      brlen = [brlen, 0];
    end;
    
    [u,f] = uniquef(anc,1);
    while (any(f>2))                                % Resolve trichotomies with zero br lengths
      i = find(f>2);                                  % Find one
      a = u(i(1));                                    % Ancestor of trichotomy
      b = brlen(a);                                   % Branch length of ancestor of trichotomy
      aa = anc(a);                                    % Ancestor of ancestor of trichotomy
      i = find(anc == a);                             % Taxa involved in trichotomy
      leni = length(i);
      i = i(leni-2:leni);                             % Reduct polychotomy to trichotomy if necessary
      next_anc = max(u)+1;                            % Next ancestor label
      anc(a) = next_anc;                              % Add new node
      brlen(a) = 0;                                   % Zero branch length leading into node
      anc(i(3)) = next_anc;
      anc = [anc, aa];
      brlen = [brlen, b];
      [u,f] = uniquef(anc,1);
    end;

    lenanc = length(anc);                           % Move root to end of ancestor function
    oldroot = find(anc==0);
    if (oldroot ~= lenanc)
      i = find(anc==oldroot);
      anc = [anc 0];
      brlen = [brlen 0];
      anc(i) = (lenanc+1)*ones(1,2);
      anc(oldroot) = [];
      i = find(anc>oldroot);
      anc(i) = anc(i)-1;
      brlen(oldroot) = [];
    end;
  end;
  
  if (~given_brlen)                               % If branch lengths weren't passed, 
    brlen = [];                                   %   don't return any
  end;

  return;
  