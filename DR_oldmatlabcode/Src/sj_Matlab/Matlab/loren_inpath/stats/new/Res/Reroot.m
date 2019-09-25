% REROOT: Given an ancestor function and vector of branch lengths, re-roots the tree 
%         with respect to a specified outgroup taxon.
%
%     [new_anc,new_brlen] = reroot(anc,brlen,outgroup,{p})
%
%           anc =       vector specifying ancestor function, with ancestor of 
%                         root specified as 0.
%           brlen =     corresponding vector of branch lengths.
%           outgroup =  index of outgroup taxon (position of outgroup in 'anc').
%           p =         optional value indicating the proportion of the distance from the 
%                         outgroup taxon to its ancestor at which the root is to be 
%                         positioned [default = 0.5].
%           ----------------------------------------------------------------------------
%           new_anc =   ancestor function of rerooted tree.
%           new_brlen = corresponding vector of branch lengths.
%

function [new_anc,new_brlen] = reroot(anc,brlen,outgroup,p)
  if (nargin < 4)
    p = [];
  end;
  if (isempty(p))
    p = 0.5;
  end;

  root = find(anc==0);                % Find current root

  [links,anc,curr_node,level] = treedivd([],anc,root,0);   % Get list of internode links
  links = links(:,1:2);               % Omit levels
  brlen = brlen(links(:,2))';         % Associate branch lengths with second col

  i = find(links(:,1)==root);         % Remove current root
  newlinks = links(i,2)';             % Add root internode to new list
  newlinks_len = sum(brlen(i));
  links(i,:) = [];
  brlen(i) = [];

  if (anc(outgroup)==root)              % If outgroup is already at the root,
    links = [links; fliplr(newlinks)];  %   restore links
    brlen = [brlen; fliplr(newlinks_len)];
    newlinks = [];
    newlinks_len = [];
  end;

  i = find(links(:,2)==outgroup);     % Create new root at outgroup
  d = links(i,1);                     % Descendant opposite outgroup

  if (isempty(i))                     % If didn't find, flip last row
    rlinks = size(links,1);
    links(rlinks,:) = fliplr(links(rlinks,:));
    i = find(links(:,2)==outgroup);     % Create new root at outgroup
    d = links(i,1);                     % Descendant opposite outgroup
  end;

  newlinks = [newlinks; outgroup root; d root];   % Add to new list
  newlinks_len = [newlinks_len; p*brlen(i); (1-p)*brlen(i)];
  links(i,:) = [];                    % Delete from old list
  brlen(i) = [];

  links = [newlinks; links];          % Concatenate remainder of links
  brlen = [newlinks_len; brlen];

  [new_anc,new_brlen,links,curr_node] = lnktoanc([],[],[links brlen],root);

  return;

