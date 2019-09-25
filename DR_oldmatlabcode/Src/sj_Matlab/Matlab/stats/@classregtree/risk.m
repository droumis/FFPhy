function r=risk(t,j)
%RISK Node risk.
%   R=RISK(T) returns an N-element vector R of the risk of the
%   nodes in the tree T, where N is the number of nodes.  The risk R(J)
%   for node J is the node error E(J) weighted by the node probability
%   P(J).  
%
%   R=RISK(T,J) takes an array J of node numbers and returns the 
%   risk values for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/ERROR, CLASSREGTREE/NODEPROB.

%   Copyright 2006-2007 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2007/01/16 03:38:07 $

if nargin>=2 && ~validatenodes(t,j)
    error('stats:classregtree:risk:InvalidNode',...
          'J must be an array of node numbers or a logical array of the proper size.');
end

if nargin<2
    r = t.risk;
else
    r = t.risk(j,:);
end
