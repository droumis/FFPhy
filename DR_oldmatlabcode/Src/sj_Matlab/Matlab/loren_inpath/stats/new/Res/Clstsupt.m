% CLSTSUPT: Finds cluster-support matrix from dendrogram topology matrix.
%
%			Usage: support = clstsupt(topology)
%
%         topology =  [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%					-----------------------------------------------------------------------
%         support =   [(n-2) x (n-1)] matrix, with one row for all but the base
%                       node, specifying group membership (support) at each node.
%

%	RE Strauss, 1/26/99

function support = clstsupt(topology)
	[r,c] = size(topology);
	n = r+1;

  support = zeros(n-2,n-1);

  for node = 1:(n-2)
    pos = 0;
    t1 = topology(node,1);
    t2 = topology(node,2);

    pos = pos+1;
    if (t1 <= n)
      support(node,pos) = t1;
    else
      mem = find(support(t1-n,:)>0);
      support(node,pos:(pos+length(mem)-1)) = support(t1-n,mem);
      pos = pos+length(mem)-1;
    end;
    
    pos = pos+1;
    if (t2 <= n)
      support(node,pos) = t2;
    else
      mem = find(support(t2-n,:)>0);
      support(node,pos:(pos+length(mem)-1)) = support(t2-n,mem);
      pos = pos+length(mem)-1;
    end;

    len = length(find(support(node,:)>0));
    support(node,1:len) = sort(support(node,1:len));
  end;

	return;
